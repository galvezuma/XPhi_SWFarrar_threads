/*
 ============================================================================
 Name        : XPhi_SWFarrar_threads.c
 Author      : SGR
 Version     :
 Copyright   : AGR-248
 Description : Smith-Waterman parallelized in C, Ansi-style
 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <fcntl.h>
#include <sys/time.h>
#include <pthread.h>
#include <immintrin.h>
#include "Splitter.h"

#define NUM_TILES 228
#define VER_TIEMPO
#define MAX_LONG_TEXTO 1024
#define MAX_LONG_PROTEINA (100*1024)

// Otras macros
#define max(x, y) (((x)>=(y))?(x):(y))
#define min(x, y) (((x)< (y))?(x):(y))

uint16_t long_sec_ref;
char * sec_ref=NULL;
uint32_t num_letras = 0;
uint8_t letras[128];
int8_t * matriz=NULL;
uint16_t long_profile;
int32_t * profile = NULL;
int coste_open_gap=10, coste_extend_gap=1, umbral=120; // See http://www.biology.wustl.edu/gcg/psiblast.html
// Declaraciones derivadas de optimizaciones
char * letras_sec_ref;
static uint8_t mejor_valor_letras[128];
static int16_t mejor_puntuacion_posible;
// Declaraciones derivadas del uso de hilos
char * filename_base_fasta;
off_t * slices; // Inicio y fin del trozo a examinar por cada hilo.

// DECLARACIÓN DE TIPOS
enum Estado {INICIO, LEYENDO_PROTEINA, FINAL};
struct Proteina {
	char nombre[MAX_LONG_TEXTO];
	uint32_t long_actual_proteina;
	char proteina[MAX_LONG_PROTEINA];
};

// DECLARACIÓN DE FUNCIONES
char * load_fasta_referencia(char* filename);
int8_t * load_matriz(char* filename);
void procesar_base_fasta(char * filename, uint16_t umbral);
static void * procesar_trozo_base_fasta(void * i);

int main(int argc, char *argv[]) {
  char nombre_fichero_matriz[MAX_LONG_TEXTO] = "BLOSUM60";
  // CONTROL DE PARÁMETROS
  int pos = 1;
  while(pos < argc - 2){
    if (!strcmp(argv[pos], "-t")) { // Leer umbral (threshold)
      umbral = atoi(argv[pos+1]);
      pos += 2;
    } else if (!strcmp(argv[pos], "-m")) { // Leer matriz
      strcpy(nombre_fichero_matriz, argv[pos+1]);
      pos += 2;
    } else if (!strcmp(argv[pos], "-g")) { // Leer costes
      coste_open_gap = atoi(argv[pos+1]);
      coste_extend_gap = atoi(argv[pos+2]);
      fprintf(stdout, "Open gap set to: %d\nExtend gap set to: %d\n", coste_open_gap, coste_extend_gap);
      pos += 3;
    } else {
        fprintf(stderr, "Invalid parameter: %s\n", argv[pos]);
        return 1;
    }
  }
  if (argc < pos + 2) {
    fprintf(stderr, "Uso: SW_threaded [-t umbral] [-m matriz] [-g coste_open_gap coste_extend_gap] referencia.fasta base.fasta\n");
    return 1;
  }
  if (((uint8_t)coste_open_gap) > 0xFF) { fprintf(stderr, "Coste de apertura de hueco (%d) demasiado grande. Su tamaño es un byte.\n", coste_open_gap); return 1; }
  if (((uint8_t)coste_extend_gap) > 0xFF) { fprintf(stderr, "Coste de extensión de hueco (%d) demasiado grande. Su tamaño es un byte.\n", coste_extend_gap); return 1; }
  if (umbral > 0xFFFF) { fprintf(stderr, "Umbral (%d) demasiado grande. El valor máximo es 0xFFFF.\n", umbral); return 1; }
  // pos apunta al nombre del fichero fasta de referencia
  sec_ref = load_fasta_referencia(argv[pos]);
  if (sec_ref == NULL) { return 1; }
  matriz = load_matriz(nombre_fichero_matriz); // Esta asignación es redundante
  if (matriz == NULL) { return 1; }
	// OPTIMIZACIONES
  	  letras_sec_ref = (char *)malloc(long_sec_ref*sizeof(char));
	// Transformación de las letras de la secuencia de referencia en posiciones de la matriz de letras.
	for(int j=0; j<long_sec_ref; j++) {
		letras_sec_ref[j] = letras[(int)(sec_ref[j])];
	}

	/* Cálculo del perfil de Farrar*/
	// LLenar el perfil sin considerar Farrar (excepto porque es más largo de lo debido).
	int segLen = (64 / sizeof(int32_t));
	long_profile = ((long_sec_ref + (63 / sizeof(int32_t))) / segLen) * segLen;
	profile = (int32_t *) _mm_malloc(num_letras*long_profile*sizeof(int32_t), 64);
	printf("Tam profile %d\n", num_letras*long_profile);
	int offset=0;
	for(int x=0; x<num_letras;x++){
		offset = x * long_profile;
		int y;
		for(y=0; y<long_sec_ref; y++){
			int8_t pos_letra_ref = letras_sec_ref[y];
			profile[offset] = matriz[pos_letra_ref*num_letras+x];
			offset++;
		}
		// Rellenar a 0 el resto
		for( ; y<long_profile; y++){
			profile[offset] = 0;
			offset++;
		}
	}

	// Barajar el perfil segun Farrar (Fig. 1 de su artículo).
	int32_t columnaAux[long_profile] __attribute__((aligned(64)));
	for(int x=0; x<num_letras; x++){
		offset = x * long_profile;
		for(int y=1; y<=long_profile; y++){
			int segLen = (64 / sizeof(int32_t));
			int numSeg = long_profile / segLen;
			int base = 1+(y-1)/segLen; // base y tt son la clave del barajar.
			int tt = (y-1) % segLen;
			columnaAux[y-1] = profile[offset + base + tt * numSeg - 1];
			//printf("orij[%d]=%d farrar[%d]=%d\n", y-1, profile[offset + y-1], base + tt * numSeg - 1, columnaAux[y-1]);
		}
		memcpy(profile+offset, columnaAux, long_profile * sizeof(int32_t));
	}
/*
	for(int x=0; x<num_letras; x++){
		offset = x * long_profile;
		printf("Letra %d\n", x);
		for(int y=0; y<long_profile; y++){
			printf("%d,  ", profile[offset + y]);
		}
		printf("\n");
	}
*/

	// BLOQUE NO USADO: Ha demostrado meter más retraso que beneficio.
	// Obtención de la mejor puntuación por letra para parar Smith-Waterman cuando se sabe que no se va a superar el umbral.
	uint8_t mejor_valor_asterisco = 0;
	// Calculamos el mejor valor por cada letra
	// Calculamos el mejor valor de la última línea de la matriz (*) para no recalcularlo múltiples veces
	for(int i=0; i<num_letras; i++)
		mejor_valor_asterisco = max(mejor_valor_asterisco, matriz[(num_letras-1)*num_letras + i]);
	// Calculamos el mejor valor por cada letra en general.
	for(int x=0; x<128; x++){
		if (letras[x] == num_letras-1) {
			mejor_valor_letras[x] = mejor_valor_asterisco;
		} else {
			uint8_t mejor_valor = 0;
			for(int i=0; i<num_letras; i++)
				mejor_valor = max(mejor_valor, matriz[letras[x]*num_letras + i]);
			mejor_valor_letras[x] = mejor_valor;
		}
	}
	// Calculamos la mejor puntuación posible en base a la secuencia de referencia y a la matriz de puntuaciones.
	mejor_puntuacion_posible = 0;
	for(int x=0; x<long_sec_ref; x++) {
		mejor_puntuacion_posible += mejor_valor_letras[(unsigned)(sec_ref[x])];
	}


	filename_base_fasta = argv[pos+1];

  procesar_base_fasta(filename_base_fasta, umbral);


  //TAREAS DE FINALIZACIÓN
	fflush(stdout);
	fflush(stderr);
	_mm_free(profile);
	free(letras_sec_ref);
	free(sec_ref);
	free(matriz);
  return 0;
}

char * load_fasta_referencia(char* filename){
	char buffer_texto[MAX_LONG_TEXTO];
	char sec_ref_local[MAX_LONG_PROTEINA];
	int num_bytes;
	FILE * fichero = fopen(filename, "r");
	if (! fichero) {
		fprintf(stderr, "%s no se puede abrir.\n", filename);
		return NULL;
	}
	char * ok = fgets(buffer_texto, MAX_LONG_TEXTO, fichero);
	if (!ok || buffer_texto[0] != '>') {
		fprintf(stderr, "No parece que %s sea un fichero FASTA.\n%s\n", filename, buffer_texto);
		fclose(fichero);
		return NULL;
	}
	// La cabecera de la secuencia de referencia no es necesario guardarla
	sec_ref_local[0] =0;
	while (fgets(buffer_texto, MAX_LONG_TEXTO, fichero)) {
		num_bytes = strlen(buffer_texto);
		// Bajo ciertas circunstancias, buffer_texto puede contener el \n o el \r. Hay que quitarlo.
		while((num_bytes > 0) && (buffer_texto[num_bytes-1] < ' ')) buffer_texto[--num_bytes] = 0;
		strcat(sec_ref_local, buffer_texto);
		long_sec_ref += num_bytes;
	}
	fclose(fichero);
	char *retorno = (char *)malloc(long_sec_ref+1);
	strcpy(retorno, sec_ref_local);
	return retorno;
}


/* Esta función utiliza en la memoria compartida:
 * un uint8_t con el número de letras del alfabeto
 * 128 uint_8 con la posición de cada letra en la matriz
 * tantos int8_t como sea (num_letras+1)*(num_letras+1)
*/
int8_t * load_matriz(char* filename) {
	int rows;
	char buffer_texto[MAX_LONG_TEXTO];

	FILE * fichero = fopen(filename, "r");
	if (!fichero) {
		fprintf(stderr, "%s no se puede abrir.\n", filename);
		return NULL ;
	}
	for (int i = 0; i < 128; i++)
		letras[i] = 255;
	rows = 0;
	while (fgets(buffer_texto, MAX_LONG_TEXTO, fichero)) {
		if (buffer_texto[0] == '#')
			continue; // Las líneas que empiezan por # son comentarios.
		if (buffer_texto[0] == ' ') { //Es el alfabeto. Se asume que es el mismo en horizontal y en vertical
			int longitud = strlen(buffer_texto);
			for (int i = 0; i < longitud; i++) {
				if (buffer_texto[i] != ' ' && buffer_texto[i] != '\n'
						&& buffer_texto[i] != '\r') { //Suponemos que es una letra
					letras[(unsigned)(buffer_texto[i])] = num_letras; // Se guarda la posición de la letra
					num_letras++;
				}
			}
			// La última letra se asume que es * (cualquier otra cosa).
			int ultima_pos = num_letras - 1;
			for (int i = 0; i < 128; i++)
				if (letras[i] == 255)
					letras[i] = ultima_pos;
			matriz = (int8_t *)malloc(num_letras*num_letras*sizeof(int8_t));
			continue;
		}
		// Si se llega aquí es porque se están leyendo las filas de la matriz
		// Los datos se graban en memoria compartida conforme se van leyendo

		//Descartamos el primer elemento (una letra en columna) y el resto vamos añadiendo
		char * linea = buffer_texto;
		char * nuevalinea = linea + 1;
		int cols = 0;
		//Mientras siga habiendo enteros (linea!=nuevalinea), añadimos a la matriz
		short fin_iteracion = 0;
		while (!fin_iteracion) {
			linea = nuevalinea;
			int8_t valor = (int8_t) strtol(linea, &nuevalinea, 0);
			if (linea != nuevalinea) {
				matriz[rows*num_letras+cols] = valor;
				cols++;
			} else {
				fin_iteracion = 1;
			}
		}
		if (cols != num_letras) {
			fprintf(stderr,
					"La fila %d de %s tiene un número de valores (%d) que no coincide con el de letras (%d).\n",
					rows, filename, cols, num_letras);
			fclose(fichero);
			free(matriz);
			return matriz=NULL;
		}
		rows++;
	}
	if (rows != num_letras) {
		fprintf(stderr,
				"El número de filas de %s es %d que no coincide con el de letras (%d).\n",
				filename, rows, num_letras);
		fclose(fichero);
		free(matriz);
		return matriz=NULL;
	}
	fclose(fichero);
	return matriz;
}




// SUBPROGRAMA DE CARGA DE SECUENCIAS DE LA BASE DE SECUENCIAS EN
// LA COLA DE LA MEMORIA COMPARTIDA.
// FINALIZA CUANDO TODAS LAS SECUENCIAS SE HAN LEIDO, SE HAN METIDO
// EN LA COLA Y SE HAN PROCESADO POR PARTE DE LOS CORES.

struct {
	uint64_t MCUPS;
} retorno_core[NUM_TILES];


void procesar_base_fasta(char * filename, uint16_t umbral){
	FILE * fichero = fopen(filename, "r");
	if (! fichero) {
		fprintf(stderr, "%s no se puede abrir.\n", filename);
		return;
	}
	fclose(fichero);

	int16_t num_threads = NUM_TILES;
	slices = split(filename, num_threads);

	pthread_t hilos[num_threads];
	for(uint64_t i=0; i<num_threads; i++)
		pthread_create(hilos+i, NULL, procesar_trozo_base_fasta, (void *)i);
	for(int i=0; i<num_threads; i++)
		pthread_join(hilos[i], NULL);
	double GCUPS = 0.0;
	for(int i=0; i<num_threads; i++)
		GCUPS += retorno_core[i].MCUPS;
	GCUPS /= 1000;
	fprintf(stdout, "GCUPS: %9.5f. ", GCUPS);
	free(slices);
}

typedef struct {
	uint64_t num_nucl_procesados;
	//
	int32_t * columnaActual_Max;
	//int32_t * columna_Up;
	int32_t * columna_Left;
	int32_t * columnaPrevia_Max;
} Objeto;

void procesar_proteina(Objeto *o, struct Proteina * p);
static inline int16_t smith_waterman_farrar(Objeto*o, char *sec_database, int16_t long_sec_database);

static void * procesar_trozo_base_fasta(void * i){
	uint32_t num_sec_leidas=0;
	enum Estado estadoActual;
	struct Proteina proteinaActual;
	char buffer_texto[MAX_LONG_TEXTO];
	uint64_t yo = (uint64_t) i;
	FILE *fp = fopen(filename_base_fasta, "r");

	fseek(fp, slices[2*yo], SEEK_SET);

	/*
	 * Preparación del estado interno del hilo
	 */
	Objeto estado_interno;
	estado_interno.num_nucl_procesados = 0;
	estado_interno.columnaActual_Max =  (int32_t *)_mm_malloc((long_profile)*sizeof(int32_t), 64);
	//estado_interno.columna_Up =   (int32_t *)_mm_malloc((long_profile+1)*sizeof(int32_t), 64);
	estado_interno.columna_Left = (int32_t *)_mm_malloc((long_profile)*sizeof(int32_t), 64);
	estado_interno.columnaPrevia_Max =  (int32_t *)_mm_malloc((long_profile)*sizeof(int32_t), 64);


	struct timeval tiempo_inicio, tiempo_final;
	gettimeofday(&tiempo_inicio, NULL);

	estadoActual = INICIO;
	while ((ftell(fp) <= slices[2*yo+1]) && fgets(buffer_texto, MAX_LONG_TEXTO, fp)) {
		int num_bytes = strlen(buffer_texto);
		// Bajo ciertas circunstancias, buffer_texto puede contener el \n o el \r. Hay que quitarlo.
		while((num_bytes > 0) && (buffer_texto[num_bytes-1] < ' ')) buffer_texto[--num_bytes] = 0;
        if (estadoActual == INICIO){
            if (buffer_texto[0] == '>'){
                strcpy(proteinaActual.nombre, buffer_texto);
				proteinaActual.long_actual_proteina = 0;
				num_sec_leidas++;
                estadoActual = LEYENDO_PROTEINA;
            } else {
                fprintf(stderr, "Formato de fichero (%s) no válido. Su contenido debe empezar por >.\n", filename_base_fasta);
                fclose(fp); return NULL;
            }
        } else { // if (estadoActual == LEYENDO_PROTEINA){
            if (buffer_texto[0] == '>'){
				if (proteinaActual.long_actual_proteina == 0) {
		            fprintf(stderr, "Error inesperado en %s. Proteina actual nula.\n", filename_base_fasta);
		            fclose(fp); return NULL;
				}
				procesar_proteina(&estado_interno, &proteinaActual);
                strcpy(proteinaActual.nombre, buffer_texto);
				proteinaActual.long_actual_proteina = 0;
				num_sec_leidas++;
            } else {
				int long_linea = strlen(buffer_texto);
                memcpy(proteinaActual.proteina + proteinaActual.long_actual_proteina, buffer_texto, long_linea);
				proteinaActual.long_actual_proteina += long_linea;
            }
        }
	}
	fclose(fp);
	// Procesar la última proteina
    if (num_sec_leidas != 0){
		if (proteinaActual.long_actual_proteina == 0) {
		    fprintf(stderr, "Error inesperado en %s. Proteina actual nula.\n", filename_base_fasta);
		    return NULL;
		}
		procesar_proteina(&estado_interno, &proteinaActual);
    } else {
        fprintf(stderr, "Fichero %s vacío.\n", filename_base_fasta);
        return NULL;
    }
    //_mm_free(estado_interno.columna_Up); estado_interno.columna_Up = NULL;
    _mm_free(estado_interno.columnaActual_Max); estado_interno.columnaActual_Max = NULL;
    _mm_free(estado_interno.columnaPrevia_Max); estado_interno.columnaPrevia_Max = NULL;
    _mm_free(estado_interno.columna_Left); estado_interno.columna_Left = NULL;


	gettimeofday(&tiempo_final, NULL);
	double tiempo_total = (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0);
	retorno_core[yo].MCUPS = ((double)(estado_interno.num_nucl_procesados*long_sec_ref))/tiempo_total/1000.0;
#ifdef VER_TIEMPO
	fprintf(stdout, "Tiempo transcurrido: %9.1fms. ", tiempo_total);
	fprintf(stdout, "GCUPS: %5.5f; ", (double)retorno_core[yo].MCUPS/1000);
#endif
	fprintf(stdout, "(%d) %d secuencias leidas.\n", yo, num_sec_leidas);
	return NULL;
}

/*
 * Parte correspondiente al procesamiento en sí.
 */

void procesar_proteina(Objeto *o, struct Proteina * p){
	uint32_t puntuacion;
	o->num_nucl_procesados += p->long_actual_proteina;
	puntuacion = smith_waterman_farrar(o, p->proteina, (int16_t)(p->long_actual_proteina));
	if (puntuacion >= umbral) {
		// Mostrar hit con nombre de secuencia.
		fprintf(stdout, "Hit. Score: %d. Secuencia: %s\n", puntuacion, p->nombre);
	}
}

inline __m512i shiftLeft(__m512i a){
	int x;
	int32_t rbuffer[64/sizeof(int32_t)] __attribute__((aligned(64)));
	// Guarda en memoria
	_mm512_store_epi32(rbuffer, a);
	// Desplaza en memoria
	for(x=0; x<64/sizeof(int32_t); x++){
		rbuffer[x] = rbuffer[x+1];
	}
	rbuffer[x] = 0;
	// Carga memoria en registro y retorna
	return _mm512_load_epi32 (rbuffer);
}

inline __m512i shiftRight(__m512i a){
	int x;
	int32_t rbuffer[64/sizeof(int32_t)] __attribute__((aligned(64)));
	// Guarda en memoria
	_mm512_store_epi32(rbuffer, a);
	// Desplaza en memoria
	for(x=64/sizeof(int32_t)-1; x>0; x--){
		rbuffer[x] = rbuffer[x-1];
	}
	rbuffer[x] = 0;
	// Carga memoria en registro y retorna
	return _mm512_load_epi32 (rbuffer);
}

inline void displayV(char * s, __m512i a){
	int x;
	int32_t rbuffer[64/sizeof(int32_t)] __attribute__((aligned(64)));
	// Guarda en memoria
	_mm512_store_epi32(rbuffer, a);
	// Desplaza en memoria
	for(x=0; x<64/sizeof(int32_t); x++){
		printf("%s[%d]=%d\n", s, x, rbuffer[x]);
	}
}

inline static int16_t smith_waterman_farrar(Objeto*o, char *sec_database, int16_t long_sec_database){
	int32_t * aux_Max;
	int16_t ret_max = 0;

	__m512i vGapOpen, vGapExtend, zero;
	__m512i vF, vH, vMax, vE_j, vAux0;
	int segLen = (64 / sizeof(int32_t));
	int numSeg = (long_sec_ref + 63 / sizeof(int32_t)) / (64 / sizeof(int32_t));

		int32_t cog[segLen] __attribute__((aligned(64)));
		int32_t ceg[segLen] __attribute__((aligned(64)));
		//
		for(int x=0;x<segLen;x++) {
			cog[x] = coste_open_gap;
			ceg[x] = coste_extend_gap;
		}
		vGapOpen = _mm512_load_epi32(cog);
		vGapExtend = _mm512_load_epi32(ceg);
		zero = _mm512_xor_epi32(zero, zero);


	vMax = _mm512_xor_epi32(vMax, vMax); // vMax = <0, 0, ..., 0>

	for(int j=0; j<long_profile; j++){
		o->columnaPrevia_Max[j] = 0;
		//o->columna_Up[j] = 0;
		o->columna_Left[j] = 0;
	}
	for(int x=0; x<long_sec_database; x++){
		// vF = <0, 0, ..., 0>
		vF = _mm512_xor_epi32(vF, vF);

		// vH = vHStore[numSeg - 1] << 1
		vH = _mm512_load_epi32(o->columnaPrevia_Max + (numSeg - 1) * segLen);
		vH = shiftRight(vH);

		//
		int8_t pos_letra_database = letras[(int)(sec_database[x])];
		//printf("Letra %d %c %d\n", x, sec_database[x], pos_letra_database);
		int32_t offset = pos_letra_database * long_profile;
		int j;
		for(j=0; j<numSeg; j++){
			// vH = vH + vProfile[letra][j]
			int32_t * valor_match = profile + offset;
			offset += segLen;
			vAux0 = _mm512_load_epi32(valor_match);
			vH = _mm512_add_epi32(vH, vAux0);

			// vMax = max(vMax, vH);
			vMax = _mm512_max_epi32(vMax, vH);

			// vE[j] = max(vH, vE[j])
			// vH = max(vH, vF)
			vE_j = _mm512_load_epi32(o->columna_Left + j*segLen);
			vH = _mm512_max_epi32(vH, vE_j);
			vH = _mm512_max_epi32(vH, vF);

			// vHStore[j] = vH
			_mm512_store_epi32(o->columnaActual_Max + j*segLen, vH);

			// vAux = vH - vGapOpen
			vAux0 = _mm512_sub_epi32(vH, vGapOpen);
			vAux0 = _mm512_max_epi32(vAux0, zero);
			// vE[j] = vE[j] - vGapExtend
			vE_j = _mm512_sub_epi32(vE_j, vGapExtend);
			vE_j = _mm512_max_epi32(vE_j, zero);
			// vE[j] = max(vE[j], vAux)
			vE_j = _mm512_max_epi32(vE_j, vAux0);
			_mm512_store_epi32(o->columna_Left + j*segLen, vE_j);
			// vF = vF - vGapExtend
			vF = _mm512_sub_epi32(vF, vGapExtend);
			vF = _mm512_max_epi32(vF, zero);
			// vF = max(vF, vAux)
			vF = _mm512_max_epi32(vF, vAux0);

			// vH = vHLoad[j]
			vH = _mm512_load_epi32(o->columnaPrevia_Max + j*segLen);
		}
		// Optimización de SWAT
		/*
		for(int x=0; x<long_profile; x++){
			printf("vMax[%d]=%d\n", x, o->columnaActual_Max[x]);
		}
		printf("Numseg: %d\n", numSeg);
		displayV("F", vF);
		*/
		// vF = vF << 1
		vF = shiftRight(vF);

		j = 0;
		do { // while(AnyElement(vF > vHStore[j] - vGapOpen
			vH = _mm512_load_epi32(o->columnaActual_Max + j*segLen);
			vAux0 = _mm512_sub_epi32(vH, vGapOpen);
			vAux0 = _mm512_max_epi32(vAux0, zero);
			__mmask16 mascara = _mm512_cmpgt_epi32_mask (vF, vAux0);
			if (mascara == 0) break;
			// vHStore[j] = max(vHStore[j], vF)
			vH = _mm512_max_epi32(vH, vF);
			_mm512_store_epi32(o->columnaActual_Max + j*segLen, vH);

			// vF = vF - vGapExtend
			vF = _mm512_sub_epi32(vF, vGapExtend);
			vF = _mm512_max_epi32(vF, zero);
			if (++j >= numSeg) {
				// vF = vF << 1
				vF = shiftRight(vF);
				j = 0;
			}

		} while(1);

		//
		aux_Max = o->columnaActual_Max;
		o->columnaActual_Max = o->columnaPrevia_Max;
		o->columnaPrevia_Max = aux_Max;
		//
	}

	int32_t max[segLen] __attribute__((aligned(64)));
	_mm512_store_epi32(max, vMax);
	for(int x=1;x<segLen;x++) {
		if(max[0] < max[x]) max[0] = max[x];
	}
	if (max[0] > 32767) max[0] = 32767;
	ret_max = max[0];
	return ret_max;
}

#ifdef SW_NORMAL
static inline int16_t smith_waterman(Objeto*o, char *sec_database, int16_t long_sec_database){
	Nodo a __attribute__((aligned(4))), b __attribute__((aligned(4)));
	Nodo *nodo_superior __attribute__((aligned(4)));
	nodo_superior = &a;
	Nodo * nodo_actual __attribute__((aligned(4)));
	nodo_actual=&b;
	Nodo * aux __attribute__((aligned(4)));
	int16_t ret_max = 0;
	for(int j=0; j<=long_sec_ref; j++){
		o->columnaActual[j].max = 0;
		o->columnaActual[j].up = coste_open_gap + coste_extend_gap;
		o->columnaActual[j].left = coste_open_gap + coste_extend_gap;
	}
	int16_t mejor_puntuacion_alcanzable = 0;
	for(int x=0; x<long_sec_database; x++) {
		mejor_puntuacion_alcanzable += mejor_valor_letras[(unsigned)(sec_database[x])];
	}
	for(int x=0; x<long_sec_database; x++){
		// Inicializar cabeza de la siguiente columna
		nodo_superior->max = 0;
		nodo_superior->up = coste_open_gap + coste_extend_gap;
		nodo_superior->left = coste_open_gap + coste_extend_gap;
		int8_t pos_letra_database = letras[(int)(sec_database[x])];
		int16_t mejor_puntuacion_vertical_posible = mejor_puntuacion_posible;
		int j;
		for(j=1; j<=long_sec_ref; j++){
			int8_t pos_letra_ref = letras_sec_ref[j-1]; // La transformación a letras[sec_ref[j-1]] se hizo justo después de cargar la secuencia;
			int8_t valor_match = matriz[pos_letra_ref*num_letras+pos_letra_database];

			int16_t result_left = o->columnaActual[j].max - o->columnaActual[j].left;
			int16_t result_up = nodo_superior->max - nodo_superior->up;

			nodo_actual->max = max(0,max(max(result_left, result_up), o->columnaActual[j-1].max+valor_match));
			nodo_actual->up = nodo_actual->max - (max(result_up, nodo_actual->max - coste_open_gap) - coste_extend_gap);
			nodo_actual->left = nodo_actual->max - (max(result_left, nodo_actual->max - coste_open_gap) - coste_extend_gap);

			ret_max = max(ret_max, nodo_actual->max);
			o->columnaActual[j-1] = *nodo_superior;
			aux = nodo_superior;
			nodo_superior = nodo_actual;
			nodo_actual = aux;
			mejor_puntuacion_vertical_posible -= mejor_valor_letras[(unsigned)(sec_ref[j-1])];
			if (ret_max + mejor_puntuacion_vertical_posible < umbral){ // Nunca se va a superar el umbral en esta columna
				j++;
				break;
			}
		}
		// Este bucle solo se ejecuta si se ha hecho un break en el anterior. Sirve para rellenar el resto de la columna
		for( ; j<=long_sec_ref; j++){
			*(uint32_t *)(o->columnaActual+j-1) = (uint32_t)0x00000000; // Porque un nodo son 4 bytes
		}
		o->columnaActual[j-1] = *nodo_superior;
		mejor_puntuacion_alcanzable -= mejor_valor_letras[(unsigned)(sec_database[x])];
		if (ret_max + mejor_puntuacion_alcanzable < umbral) // Nunca se va a superar el umbral
			break;
	}
	return ret_max;
}
#endif

#ifdef SW_FULL
// Versión con dos columnas. Se guardan los datos completos y se usa
// optimización SWAT y de profile.
inline static int16_t smith_waterman(Objeto*o, char *sec_database, int16_t long_sec_database){
	Nodo * aux __attribute__((aligned(8)));
	int16_t ret_max = 0;
	for(int j=0; j<=long_sec_ref; j++){
		o->columnaPrevia[j].max = 0;
		o->columnaPrevia[j].up = 0;
		o->columnaPrevia[j].left = 0;
	}
	for(int x=0; x<long_sec_database; x++){
		// Inicializar cabeza de la siguiente columna
		o->columnaActual[0].max = 0;
		o->columnaActual[0].up = 0;
		o->columnaActual[0].left = 0;
		int8_t pos_letra_database = letras[(int)(sec_database[x])];
		int16_t offset = pos_letra_database * long_profile;
		int j;
		for(j=1; j<=long_sec_ref; j++){
			// Sustituido por optimización Profile
			//int8_t pos_letra_ref = letras_sec_ref[j-1]; // La transformación a letras[sec_ref[j-1]] se hizo justo después de cargar la secuencia;
			//int8_t valor_match = matriz[pos_letra_ref*num_letras+pos_letra_database];
			int8_t valor_match = profile[offset];
			offset++;

			// Sustituido por Optimización de SWAT
			// o->columnaActual[j].up = max(o->columnaActual[j-1].up, o->columnaActual[j-1].max - coste_open_gap) - coste_extend_gap;

			////o->columnaActual[j].left = max(o->columnaPrevia[j].left, o->columnaPrevia[j].max - coste_open_gap) - coste_extend_gap;

			o->columnaActual[j].max = max(0, o->columnaPrevia[j-1].max+valor_match);
			////o->columnaActual[j].max = max(0, max(o->columnaPrevia[j-1].max+valor_match, o->columnaActual[j].left));
			// Sustituido por Optimización de SWAT
			// o->columnaActual[j].max = max(0,max(max(o->columnaActual[j].up, o->columnaActual[j].left), o->columnaPrevia[j-1].max+valor_match));

			ret_max = max(ret_max, o->columnaActual[j].max);
		}
		// Optimización de SWAT
		for(j=1; j<=long_sec_ref; j++){
			if (o->columnaPrevia[j].left != 0 || o->columnaPrevia[j].max > coste_open_gap + coste_extend_gap){
				o->columnaActual[j].left = max(o->columnaPrevia[j].left, o->columnaPrevia[j].max - coste_open_gap) - coste_extend_gap;
				o->columnaActual[j].max = max(o->columnaActual[j].max, o->columnaActual[j].left);
			} else {
				o->columnaActual[j].left = 0;
			}
			if (o->columnaActual[j-1].up != 0 || o->columnaActual[j-1].max > coste_open_gap + coste_extend_gap){
				o->columnaActual[j].up = max(o->columnaActual[j-1].up, o->columnaActual[j-1].max - coste_open_gap) - coste_extend_gap;
				o->columnaActual[j].max = max(o->columnaActual[j].max, o->columnaActual[j].up);
			} else o->columnaActual[j].up = 0;
		}
		aux = o->columnaActual;
		o->columnaActual = o->columnaPrevia;
		o->columnaPrevia = aux;
	}
	return ret_max;
}
#endif
