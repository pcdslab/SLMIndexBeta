// Copyright 2012 Eric Liang

#ifndef BWT_H_
#define BWT_H_

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

typedef unsigned char Peak;

// Structure to store data of a rotation
struct rotation
{
    unsigned int index;
    Peak *suffix;

    struct rotation operator=(struct rotation a)
    {
        index=a.index;
        suffix=a.suffix;
        return *this;
    }
};

unsigned int *computeSuffixArray(Peak *input_text, unsigned int len_text);
unsigned char *findLastChar(Peak *input_text, unsigned int *suffix_arr, unsigned int n);
unsigned char *findFirstChar(Peak *input_text, unsigned int *suffix_arr, unsigned int n);
void BWTmergeSort(struct rotation arr[], int left, int right);
void BWTmerge(struct rotation arr[], int left, int middle, int right);

#endif  // BWT_H_
