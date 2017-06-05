#include <stdio.h>
#include <cstring> // for strcmp
#include "assert.h"
#include "dog_str.h"

#define NDIMS 0  // Why is this here? (-DS)

#include "tensors1d.h"
#include "debug.h"

bool str_eq(const char* str1, const char* str2)
{
    // for case-insensitive matching change this to strcasecmp
    return !strcmp(str1,str2);
}

// convert comma-separated list to array of integers.
//
// user is responsible to
// delete [] int_array;
//
int new_array_from_str(
        const char* str_in,
        int* &int_array,
        int first_idx,
        char field_sep)
{
    assert(first_idx==0 || first_idx==1);

    int_array = 0;
    char *str = (char*)str_in;
    char *ptr = str;
    if(ptr==0) return -1;
    if(*ptr==0) return 0;

    // allocate int_array
    int number_of_numbers=0;
    {
        // count number of commas in string
        int number_of_commas=0;
        for(char *ptr=str; *ptr; ptr++) if(*ptr == field_sep) number_of_commas++;
        number_of_numbers = number_of_commas+1;
        int array_size = first_idx+number_of_numbers;
        int_array=new int[array_size];
        for(int i=0;i<array_size;i++) int_array[i]=0;
    }

    // populate int_array with numbers from list
    ptr=str;
    int idx=first_idx;
    while(*ptr) 
    {
        int val;
        int successful_match = sscanf(ptr,"%d",&val);
        if(!successful_match) {
            printf("problem: could not match %s as integer\n",ptr);
            delete [] int_array; int_array=0; return -1;
        } else {
            int_array[idx++]=val;
        }
        // advance ptr past next comma
        while(*ptr && *ptr++ != field_sep);
    }
    assert(idx==(first_idx+number_of_numbers));
    return number_of_numbers;
}

// int count_fields(const char* str_in, char field_sep)
// {
//   char *str = (char*)str_in;
//   char *ptr = str;
//   // a NULL or empty string is considered to have no fields
//   if(ptr==0) return 0;
//   if(*ptr==0) return 0;
//   int num_fields=1;
//   // count number of field_sep in string
//   for(char *ptr=str; *ptr; ptr++) if(*ptr == field_sep) num_fields++;
//   return num_fields;
// }
// 
// int count_char(const char* str, char c)
// {
//   if(str==0) return -1;
//   int num_c=0;
//   // count number of occurrences of c in string
//   for(char *ptr=(char*)str; *ptr; ptr++) if(*ptr == c) num_c++;
//   return num_c;
// }

// field_sep is typically ',' or '\n'
bool str_into_tensor(const char* str, iTensorBase& t, char field_sep)
{
    int idx=0;
    bool ret=true;
    char* ptr = (char*)str;
    while(*ptr && idx < t.numel()) 
    {
        int val;
        int successful_match = sscanf(ptr,"%d",&val);
        if(!successful_match) {
            Wprintf("could not match %s as integer\n",ptr);
            ret=false;
        } else {
            t.vset(idx++, val);
        }
        // advance ptr past next comma
        while(*ptr && *ptr++ != field_sep);
    }
    return ret;
};

// field_sep is typically ',' or '\n'
bool str_into_tensor(const char* str, dTensorBase& t, char field_sep)
{
    int idx=0;
    bool ret=true;
    char* ptr = (char*)str;
    while(*ptr && idx < t.numel()) 
    {
        double val;
        int successful_match = sscanf(ptr,"%lf",&val);
        if(!successful_match) {
            Wprintf("could not match %s as double\n",ptr);
            ret=false;
        } else {
            t.vset(idx++, val);
        }
        // advance ptr past next comma
        while(*ptr && *ptr++ != field_sep);
    }
    return ret;
};

