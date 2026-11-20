#ifndef JULIA_TEST_PKG_H
#define JULIA_TEST_PKG_H

#ifdef __cplusplus
extern "C" {
#endif

struct MyStruct {
    int a;
    double b;
};

void printint(int n);
float mean_of_two(int x, int y);

int struct_test(MyStruct s);

#ifdef __cplusplus
}
#endif

#endif
