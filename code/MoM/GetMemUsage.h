#ifndef GETMEMUSAGE_H
#define GETMEMUSAGE_H
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
    
// written by Maxime Martinasso, november 2011

static long parseLine(char* line)
{
    long i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atol(line);
    return i;
}
    

long MemoryUsageGetPeak()
{ 
    FILE* file = fopen("/proc/self/status", "r");
    long result = 0;
    char line[128];
    
    while (fgets(line, 128, file) != NULL)
    {
        if (strncmp(line, "VmPeak:", 6) == 0) 
        {
           result = parseLine(line);
           break;
        }
    }
    fclose(file);
    return result*1024;
}

#endif
