#include <stdio.h>
#include <stdlib.h>
#include "static_array.h"

TYPE *static_create(int dim){
    TYPE *arr;
    arr = calloc(dim, sizeof(TYPE));

    return arr;
}

void static_random_filler(TYPE *v, int dim){
    int i;
    for(i = 0; i < dim; i++){
        v[i] = rand()%50;
    }
}

/* int v_get_size(TYPE *v){
    int dim;
    dim = sizeof(v) / sizeof(v[0]);

    return dim;
}*/

void static_print_union(TYPE arr1[], TYPE arr2[], int len1, int len2) {

    int f/* frequenza */, i /* primo array */, j/* secondo array */, k = 0/* terzo array */;
    TYPE arr3[len1 + len2]/* array dove memorizzare l'unione */;

    /* copio il primo array nel nuovo array */
    for (i = 0; i < len1; i++) {
        arr3[k] = arr1[i];
        k++;
    }

    /* ciclo il secondo array */
    for (i = 0; i < len2; i++) {
        /* frequenza valore */
        f = 0;
        /* ciclo il primo array e lo comparo al secondo */
        for (j = 0; j < len1; j++) {
            /* se il valore è presente in entrambi aumento la freuenza del contatore */
            if (arr2[i] == arr1[j]) {
                f = 1;/* vuol dire che il valore è comparso in entrambi gli array */
            }
        }
        /* se f = 0 vuol dire che il valore non c'è nel primo array, ma solo nel secondo */
        if (f == 0) {
            /* copio il valore dal secondo array al terzo e incremento k */
            arr3[k] = arr2[i];
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

void static_print_intersection(TYPE *arr1, TYPE *arr2, int len1, int len2) {
    TYPE arr3[len1+len2];
    int i/* primo array */, j/*secondo array*/, k = 0/* terzo array */;

    /* ciclo il primo array */
    for (i = 0; i < len1; i++) {
        /* ciclo il secondo array */
        for (j = 0; j < len2; j++) {
            /* se il valore del secondo array è presente nel primo allora lo copia in arr3 */
            if (arr1[i] == arr2[j]) {
                arr3[k] = arr1[i];
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
    free(arr3);
}

void static_del_at_pos(TYPE *v, int dim, int pos){
    if(pos < 0 || pos > dim)
        printf("Posizione inesistente!");
    else
        v[pos-1] = 0;
}

void modify_at_pos(TYPE *v, int dim, int pos, TYPE value){
    if(pos < 0 || pos > dim)
        printf("Posizione inesistente!");
    else
        v[pos] = value;
}

void static_del_value(TYPE *v, int dim, TYPE value){
    int i;

    for(i = 0; i < dim; i++){
        if(v[i] == value)
            break;
    }
    v[i] = 0;
}

void static_modify_value(TYPE *v, int dim, TYPE previous, TYPE value){
    int i;

    for(i = 0; i < dim; i++){
        if(v[i] == previous)
            v[i] = value;
    }
}

TYPE *static_invert_v(TYPE *v, int dim){
    int i, j = 0;
    TYPE *array;
    array = static_create(dim);

    for (i = dim-1; i >= 0; i--, j++){
        array[j] = v[i];
    }

    return array;
}

void static_del(TYPE *v){
    free(v);
}

int static_get_value_pos(TYPE *v, int dim, TYPE value){
    int i;
    for(i = 0; i < dim; i++){
        if(v[i] == value)
            break;
    }
    return i;
}

void static_print(TYPE *v, int dim){
    int i;

    printf("\n Array di %d elementi\n", dim);
    printf("+----------------------------+\n");
    for(i = 0; i < dim; i++){
        printf(" Valore in posizione %d: %d\n", i, v[i]);
    }
}



/*int main() {
    int *arr1, *arr2;
    arr1 = static_create(10);
    //arr2 = static_create(5);

    static_random_filler(arr1, 10);
    //static_random_filler(arr2, 5);

    modify_at_pos(arr1, 10, 3, 99);
    v_print(arr1, 10);
    del_value(arr1, 10, 99);
    v_print(arr1, 10);

    printf("\n\nValue pos: %d", get_value_pos(arr1, 10, 28));

    static_del(arr1);
    return 0;
}*/