#ifndef SIM_BP_H
#define SIM_BP_H
#include <math.h>

typedef struct bp_params{
    unsigned long int K;
    unsigned long int M1;
    unsigned long int M2;
    unsigned long int N;
    char*             bp_name;
}bp_params;

// Put additional data structures here as per your requirement

class counter_array
{
    public:
        int total_rows;
        int* counter;

        counter_array(int M2)
        {
            total_rows = (int)(pow(2,M2));
            counter = new int[total_rows]();

            for(int i=0;i<total_rows;i++)
            {
                counter[i] = 2;
            }
        }
};


class chooser_array
{
    public:
        int total_rows;
        int* chooser;

        chooser_array(int K)
        {
            total_rows = (int)(pow(2,K));
            chooser = new int[total_rows]();

            for(int i=0;i<total_rows;i++)
            {
                chooser[i] = 1;
            }
        }
};

typedef struct global_register
{
    unsigned long int reg = 0;

}global_reg;

#endif
