#ifndef __DOG_STR_H__
#define __DOG_STR_H__

class iTensorBase;

// test strings for equality
bool str_eq(const char* str1, const char* str2);

// copy the contents of input string into the iTensor
bool str_into_tensor(const char* str, iTensorBase& t, char sep=',');

// Create a new array redirect the input pointer to this array from input string
// User is expected to delete the allocated memory when finished
int new_array_from_str( const char* str, int* &arr, int first_idx, char sep=',');

#endif
