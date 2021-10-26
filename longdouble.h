#include <stdio.h>

typedef long double longdouble;

//int snprintf(char *restrict str, size_t size,
//                   const char *restrict format, ...);

int main(void) {
    long double as = 123;

    int showpoint = 1;
    int prec = 32;
    int left = 1;
    int width = 50;

    printf("%#-*.*Lf\n", width, prec, as);

    return 0;
}
