#include <iostream>
#include "static_array.h"

TYPE *static_create(int dim){
    TYPE *arr;
    arr = new TYPE[dim];// (dim, sizeof(TYPE));

    return arr;
}

TYPE **static_pointer_create(int dim)
{
    TYPE **arr;
    arr = new TYPE*[dim]; // (dim, sizeof(TYPE));

    return arr;
}

TYPE* static_matrix_create(int d = 10, int r = 10, int c = 10, float v = 0.0){
    float **data= static_pointer_create(d);

    for(int i = 0; i < d; i++){
        data[i] = static_create(r*c);
    }

    return *data;
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
    delete v;
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
    float *vect;

    vect = static_create(10);
    static_print(vect, 10);
    stat
    return 0;
}*/