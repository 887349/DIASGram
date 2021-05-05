#include <stdio.h>
#include <stdlib.h>
#include "dynamic_array.h"



dynamic_array_t *dynamic_create(int dim){
    dynamic_array_t *v;
    v = (dynamic_array_t*)malloc(sizeof(dynamic_array_t));
    v->values = (TYPE*)calloc(dim, sizeof(TYPE));
    v->available = (unsigned int)dim;
    v->used = 0;

    return v;
}

void dynamic_random_filler(dynamic_array_t *v){
    int i;
    for(i = 0; i < v->available; i++){
        v->values[i] = rand()%50;
    }
    v->used = v->available;
}

void dynamic_print(dynamic_array_t *v){
    unsigned int dim, i;
    dim = v->available;
    printf("\n Array di %d elementi\n", dim);
    printf("+------------------------------+\n");
    for(i = 0; i < dim; i++){
        printf(" Valore in posizione %d: %d\n", i, v->values[i]);
    }
}

void dynamic_del(dynamic_array_t *v){
    free(v->values);
    free(v);
}

void dynamic_sharp_amp(dynamic_array_t *v){
    unsigned int i;
    TYPE *arr;
    arr = (TYPE*)calloc(v->available + 1, sizeof(TYPE));

    for(i = 0; i < v->available; i++){
        arr[i] = v->values[i];
    }
    free(v->values);
    v->values = arr;
    v->available++;
}

void dynamic_amp(dynamic_array_t *v){
    unsigned int i;
    TYPE *arr;
    arr = (TYPE*)calloc(v->available + 5, sizeof(TYPE));

    for(i = 0; i < v->available; i++){
        arr[i] = v->values[i];
    }
    free(v->values);
    v->values = arr;
    v->available += 5;
}

void dynamic_bulk_amp(dynamic_array_t *v){
    unsigned int i;
    TYPE *arr;
    arr = (TYPE*)calloc(v->available*2, sizeof(TYPE));

    for(i = 0; i < v->available; i++){
        arr[i] = v->values[i];
    }
    free(v->values);
    v->values = arr;
    v->available *= 2;
}

void dynamic_add(dynamic_array_t *v, TYPE value){
    unsigned int i;
    if(v->used == v->available)
        dynamic_sharp_amp(v);
    i = (v->used);
    v->values[i] = value;
    v->used++;
}

void dynamic_modify_value(dynamic_array_t *v, int previous, int value){
    int i;
    for(i = 0; i < v->used; i++){
        if(v->values[i] == previous)
            v->values[i] = value;
    }
}

void dynamic_modify_at_position(dynamic_array_t *v, unsigned int pos, TYPE value){
    if(pos < 0 || pos > v->available)
        printf("\nPosizione inesistente!\n");
    v->values[pos] = value;
}

unsigned int dynamic_get_position(dynamic_array_t *v, TYPE value){
    unsigned int i;
    for(i = 0; i < v->used; i++){
        if(v->values[i] == value)
            break;
    }
    return i;
}

void dynamic_print_intersection(dynamic_array_t *arr1, dynamic_array_t *arr2){
    unsigned int len1 = arr1->used, len2 = arr2->used;
    TYPE arr3[len1+len2];
    int i/* primo array */, j/*secondo array*/, k = 0/* terzo array */;

    /* ciclo il primo array */
    for (i = 0; i < len1; i++) {
        /* ciclo il secondo array */
        for (j = 0; j < len2; j++) {
            /* se il valore del secondo array è presente nel primo allora lo copia in arr3 */
            if (arr1->values[i] == arr2->values[j]) {
                arr3[k] = arr1->values[i];
                /* incremento k manualmente solo se viene copiato un valore */
                k++;
            }
        }
    }


    printf("\nIntersection of two array is:");
    if(k == 0) {
        printf("none\n");
    }else {
        for (i = 0; i < k; i++) {
            printf("%d ", arr3[i]);
        }
    }
}

void dynamic_print_union(dynamic_array_t *arr1, dynamic_array_t *arr2){
    unsigned int len1 = arr1->used, len2 = arr2->used;
    int f/* frequenza */, i /* primo array */, j/* secondo array */, k = 0/* terzo array */;
    TYPE arr3[len1 + len2]/* array dove memorizzare l'unione */;

    /* copio il primo array nel nuovo array */
    for (i = 0; i < len1; i++) {
        arr3[k] = arr1->values[i];
        k++;
    }

    /* ciclo il secondo array */
    for (i = 0; i < len2; i++) {
        /* frequenza valore */
        f = 0;
        /* ciclo il primo array e lo comparo al secondo */
        for (j = 0; j < len1; j++) {
            /* se il valore è presente in entrambi aumento la freuenza del contatore */
            if (arr2->values[i] == arr1->values[j]) {
                f = 1;/* vuol dire che il valore è comparso in entrambi gli array */
            }
        }
        /* se f = 0 vuol dire che il valore non c'è nel primo array, ma solo nel secondo */
        if (f == 0) {
            /* copio il valore dal secondo array al terzo e incremento k */
            arr3[k] = arr2->values[i];
            k++;
        }
        /* ricomincio il ciclo e prendo in considerazione il valore dop odel secondo array
         * e lo comparo con tutti i valori del primo */
    }

    /* stampo l'array risultante */
    printf("\nUnion of two array is: ");
    for (i = 0; i < k; i++) {
        printf("%d ", arr3[i]);
    }
}

void dynamic_del_value(dynamic_array_t *v, TYPE value){
    int i;
    for(i = 0; i < v->used; i++) {
        if (v->values[i] == value)
            break;
    }
    v->values[i] = 0;
    v->used--;
}

/*int main() {
    dynamic_array_t *arr1, *arr2;
    arr1 = dynamic_create(10);
    arr2 = dynamic_create(5);

    dynamic_random_filler(arr1);
    dynamic_random_filler(arr2);

    dynamic_modify_at_position(arr1, 3, 100);
    dynamic_modify_at_position(arr2, 3, 100);

    dynamic_modify_at_position(arr1, 2, 555);
    dynamic_modify_at_position(arr2, 4, 555);



    dynamic_print(arr1);
    dynamic_print(arr2);

    dynamic_print_union(arr1, arr2);



    dynamic_del(arr1);
    dynamic_del(arr2);
    return 0;
}*/