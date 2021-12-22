#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sim_bp.h"


/*  argc holds the number of command line arguments
    argv[] holds the commands themselves

    Example:-
    sim bimodal 6 gcc_trace.txt
    argc = 4
    argv[0] = "sim"
    argv[1] = "bimodal"
    argv[2] = "6"
    ... and so on
*/
unsigned int number_of_predictions = 0;
unsigned int number_of_mispredictions = 0;

float branch_misprediction_rate = 0.0;
bool verbose = false;

void bimodal_branch_predictor(bp_params parameters,counter_array* counter_arr,unsigned long int addr,char outcome)
{
    
    int m = parameters.M2;
    number_of_predictions += 1;

    unsigned int PC = ((addr >> 2) & (counter_arr->total_rows - 1));

    unsigned int prediction = counter_arr->counter[PC];
    bool taken = false;
    // 0 -> not taken, 1-> taken
    if(prediction >= 2)
        taken = true;

    // true_outcome_taken -> 0 not taken(original outcome) ,  true_outcome_taken -> 1 taken(original outcome)
    bool true_outcome_taken = false;
    if(outcome == 't')
        true_outcome_taken = true;

    if(taken != true_outcome_taken)
    {
        number_of_mispredictions += 1;
    }

    int counter_val = counter_arr->counter[PC];

    if(verbose)
        printf("\n BP: %d   %d\n",PC,counter_val);

    if(true_outcome_taken)
    {
        //branch is taken actually
        //increment counter entry
        if(counter_arr->counter[PC] < 3)
        {
            counter_arr->counter[PC] += 1;
        }
    }

    else
    {
        //branch is not taken actually
        //decrement counter entry
        if(counter_arr->counter[PC] > 0)
        {
            counter_arr->counter[PC] -= 1;
        }
    }
    counter_val = counter_arr->counter[PC];

    if(verbose)
        printf("\n BU: %d   %d\n",PC,counter_val);

    return;

}


void gshare_branch_predictor(bp_params parameters,global_reg *reg,counter_array* counter_arr,unsigned long int addr,char outcome)
{
    //printf("\n HERE1 \n");
    int m1 = parameters.M1;
    int n = parameters.N;

    int n_bits = ((int)(pow(2,n)));
    int global_history_index = ((reg->reg) & (n_bits - 1));
    //printf("\nGLOB_%d\n",reg->reg);
    unsigned int PC = ((addr >> 2) & (counter_arr->total_rows - 1));
    int mn_diff = m1 - n;
    int mn_pow = ((int)(pow(2,mn_diff)));
    int lower_bits_addr = (PC) & (mn_pow - 1);
    //printf("\n LOW_%d\n",lower_bits_addr);
    int upper_bits_addr = (PC >> mn_diff);
    //printf("\n UPP_%d\n",upper_bits_addr);
    int xored = global_history_index ^ upper_bits_addr;

    int PC_index = (xored << mn_diff) | lower_bits_addr;
   // printf("\nCOUNT_%d\n",PC_index);
    unsigned int prediction = counter_arr->counter[PC_index];
    //printf("\n HERE3 \n");
    bool taken = false;

    number_of_predictions += 1;

    if(verbose)
        printf("\n Global History reg:   %d",global_history_index);

    if(prediction >= 2)
        taken = true;

    bool true_outcome_taken = false;

    if(outcome == 't')
        true_outcome_taken = true;

    if(taken != true_outcome_taken)
    {
        number_of_mispredictions += 1;
    }
    
    int counter_val = counter_arr->counter[PC_index];
   // printf("\n HERE\n");
    if(verbose)
        printf("\n GP: %d   %d\n",PC_index,counter_val);

    if(true_outcome_taken)
    {
        if(counter_arr->counter[PC_index] < 3)
        {
            counter_arr->counter[PC_index] += 1;
        }
        //Update the global history register
        reg->reg  = (reg->reg >> 1) | (1 << (n-1));
        //printf("\nREG_%d\n",reg->reg);
    }

    else
    {
        if(counter_arr->counter[PC_index] > 0)
        {
            counter_arr->counter[PC_index] -= 1;
        }

        reg->reg  = (reg->reg >> 1) | (0 << (n-1));
    }

    counter_val = counter_arr->counter[PC_index];

    if(verbose)
    {
        printf("\n GU: %d   %d\n",PC_index,counter_val);
    }
    return;
}


void hybrid_branch_predictor(bp_params parameters,global_reg *reg,counter_array* counter_arr_bimodal,counter_array* counter_arr_gshare,chooser_array *chooser_arr,unsigned long int addr,char outcome)
{
    //choosing step
    int k = parameters.K;
    //printf("\nK_%d",k);
    int n = parameters.N;
    number_of_predictions += 1;
    unsigned int choose_PC = ((addr >> 2) & (chooser_arr->total_rows - 1));
    unsigned int choose_prediction = chooser_arr->chooser[choose_PC];
    bool true_outcome_taken = false;

    if(outcome == 't')
        true_outcome_taken = true;

    //gshare predict
    int m1 = parameters.M1;
        

    int n_bits = ((int)(pow(2,n)));
    int global_history_index = ((reg->reg) & (n_bits - 1));
    //printf("\nGLOB_%d\n",reg->reg);
    unsigned int PC = ((addr >> 2) & (counter_arr_gshare->total_rows - 1));
    int mn_diff = m1 - n;
    int mn_pow = ((int)(pow(2,mn_diff)));
    int lower_bits_addr = (PC) & (mn_pow - 1);
    //printf("\n LOW_%d\n",lower_bits_addr);
    int upper_bits_addr = (PC >> mn_diff);
    //printf("\n UPP_%d\n",upper_bits_addr);
    int xored = global_history_index ^ upper_bits_addr;

    int PC_index = (xored << mn_diff) | lower_bits_addr;
   // printf("\nCOUNT_%d\n",PC_index);
    unsigned int prediction = counter_arr_gshare->counter[PC_index];
    //printf("\n HERE3 \n");
    bool taken_gshare = false;

        
    if(verbose)
        printf("\n Global History reg:   %d",global_history_index);

    if(prediction >= 2)
        taken_gshare = true;    


    //bimodal predict   
    int m = parameters.M2;
    unsigned int PC_bimodal = ((addr >> 2) & (counter_arr_bimodal->total_rows - 1));
    unsigned int prediction_bimodal = counter_arr_bimodal->counter[PC_bimodal];
    bool taken_bimodal = false;
    if(prediction_bimodal >= 2)
        taken_bimodal = true;  

    if(choose_prediction >= 2)
    {
        //gshare update
        if(taken_gshare != true_outcome_taken)
        {
            number_of_mispredictions += 1;
        }
        
        int counter_val = counter_arr_gshare->counter[PC_index];
    // printf("\n HERE\n");
        if(verbose)
            printf("\n GP: %d   %d\n",PC_index,counter_val);

        if(true_outcome_taken)
        {
            if(counter_arr_gshare->counter[PC_index] < 3)
            {
                counter_arr_gshare->counter[PC_index] += 1;
            }
            
        }

        else
        {
            if(counter_arr_gshare->counter[PC_index] > 0)
            {
                counter_arr_gshare->counter[PC_index] -= 1;
            }
        
        }

        counter_val = counter_arr_gshare->counter[PC_index];

        if(verbose)
        {
            printf("\n GU: %d   %d\n",PC_index,counter_val);
        }
        
        
    }
    else
    {
        //bimodal update
        if(taken_bimodal != true_outcome_taken)
        {
            number_of_mispredictions += 1;
        }

        int counter_val = counter_arr_bimodal->counter[PC_bimodal];

        if(verbose)
            printf("\n BP: %d   %d\n",PC,counter_val);

        if(true_outcome_taken)
        {
            //branch is taken actually
            //increment counter entry
            if(counter_arr_bimodal->counter[PC_bimodal] < 3)
            {
                counter_arr_bimodal->counter[PC_bimodal] += 1;
            }
        }

        else
        {
            //branch is not taken actually
            //decrement counter entry
            if(counter_arr_bimodal->counter[PC_bimodal] > 0)
            {
                counter_arr_bimodal->counter[PC_bimodal] -= 1;
            }
        }
        counter_val = counter_arr_bimodal->counter[PC_bimodal];

        if(verbose)
            printf("\n BU: %d   %d\n",PC_bimodal,counter_val);

    }

    if(true_outcome_taken)
    {
        //Update the global history register
        reg->reg  = (reg->reg >> 1) | (1 << (n-1));
        //printf("\nREG_%d\n",reg->reg);
    }

    else
    {
        reg->reg  = (reg->reg >> 1) | (0 << (n-1));
    }
    //Updating chooser table

    if((taken_bimodal!=true_outcome_taken)&&(taken_gshare==true_outcome_taken))
    {
        //gshare correct,bimodal incorrect
        if(chooser_arr->chooser[choose_PC] < 3)
        {
            chooser_arr->chooser[choose_PC] += 1;
        }

    }

    else if((taken_bimodal==true_outcome_taken)&&(taken_gshare!=true_outcome_taken))
    {
        //bimodal correct,gshare incorrect
        if(chooser_arr->chooser[choose_PC] > 0)
        {
            chooser_arr->chooser[choose_PC] -= 1;
        }

    }

    else
    {
        //no action
    }

    return;


}

int main (int argc, char* argv[])
{
    FILE *FP;               // File handler
    char *trace_file;       // Variable that holds trace file name;
    bp_params params;       // look at sim_bp.h header file for the the definition of struct bp_params
    char outcome;           // Variable holds branch outcome
    unsigned long int addr; // Variable holds the address read from input file
    
    counter_array *counter;  //gshare
    counter_array *counter1; //bimodal
    chooser_array *chooser;

    if (!(argc == 4 || argc == 5 || argc == 7))
    {
        printf("Error: Wrong number of inputs:%d\n", argc-1);
        exit(EXIT_FAILURE);
    }
    
    params.bp_name  = argv[1];
    
    // strtoul() converts char* to unsigned long. It is included in <stdlib.h>
    if(strcmp(params.bp_name, "bimodal") == 0)              // Bimodal
    {
        if(argc != 4)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M2       = strtoul(argv[2], NULL, 10);
        trace_file      = argv[3];
        printf("COMMAND\n%s %s %lu %s\n", argv[0], params.bp_name, params.M2, trace_file);

        counter = new counter_array(params.M2);
    }
    else if(strcmp(params.bp_name, "gshare") == 0)          // Gshare
    {
        if(argc != 5)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M1       = strtoul(argv[2], NULL, 10);
        params.N        = strtoul(argv[3], NULL, 10);
        trace_file      = argv[4];
        printf("COMMAND\n%s %s %lu %lu %s\n", argv[0], params.bp_name, params.M1, params.N, trace_file);
        counter = new counter_array(params.M1);

    }
    else if(strcmp(params.bp_name, "hybrid") == 0)          // Hybrid
    {
        if(argc != 7)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.K        = strtoul(argv[2], NULL, 10);
        params.M1       = strtoul(argv[3], NULL, 10);
        params.N        = strtoul(argv[4], NULL, 10);
        params.M2       = strtoul(argv[5], NULL, 10);
        trace_file      = argv[6];
        printf("COMMAND\n%s %s %lu %lu %lu %lu %s\n", argv[0], params.bp_name, params.K, params.M1, params.N, params.M2, trace_file);
        counter = new counter_array(params.M1);
        counter1 = new counter_array(params.M2);
        chooser = new chooser_array(params.K);

    }
    else
    {
        printf("Error: Wrong branch predictor name:%s\n", params.bp_name);
        exit(EXIT_FAILURE);
    }
    
    // Open trace_file in read mode
    FP = fopen(trace_file, "r");
    if(FP == NULL)
    {
        // Throw error and exit if fopen() failed
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }
    
    char str[2];
    int c = 0;
    printf("OUTPUT\n");
    global_reg register1;
    while(fscanf(FP, "%lx %s", &addr, str) != EOF)
    {
        
        outcome = str[0];
        if(verbose)
        {
            printf("\n =%d\t",c++);
            if (outcome == 't')
                printf("%lx %s\n", addr, "t");           // Print and test if file is read correctly
            else if (outcome == 'n')
                printf("%lx %s\n", addr, "n");          // Print and test if file is read correctly
        }
        
        /*************************************
            Add branch predictor code here
        **************************************/
       if(strcmp(params.bp_name, "bimodal") == 0)
       {
         bimodal_branch_predictor(params,counter,addr,outcome);  
       }
       else if(strcmp(params.bp_name, "gshare") == 0) 
       {
           gshare_branch_predictor(params,&register1,counter,addr,outcome);
       } 
       else if(strcmp(params.bp_name, "hybrid") == 0)
       {
           hybrid_branch_predictor(params,&register1,counter1,counter,chooser,addr,outcome);
       }
       
    }
    printf("number of predictions:   %d \n",number_of_predictions);
    printf("number of mispredictions: %d  \n",number_of_mispredictions);
    branch_misprediction_rate = (((float)(number_of_mispredictions))/((float)(number_of_predictions))) * 100;
    printf("misprediction rate:       %.2f%c\n",branch_misprediction_rate,'%');

    if(strcmp(params.bp_name, "bimodal") == 0)
    {
        printf("FINAL BIMODAL CONTENTS");
        for(int i=0;i<counter->total_rows;i++)
        {
            printf("\n%d  %d",i,counter->counter[i]);
        }
    }
    else if(strcmp(params.bp_name, "gshare") == 0) 
    {
        printf("FINAL GSHARE CONTENTS");
        for(int i=0;i<counter->total_rows;i++)
        {
            printf("\n%d  %d",i,counter->counter[i]);
        }
    } 
    else if(strcmp(params.bp_name, "hybrid") == 0)
    {
        printf("FINAL CHOOSER CONTENTS");
        for(int i=0;i<chooser->total_rows;i++)
        {
            printf("\n%d  %d",i,chooser->chooser[i]);
        }
        printf("\nFINAL GSHARE CONTENTS");
        for(int i=0;i<counter->total_rows;i++)
        {
            printf("\n%d  %d",i,counter->counter[i]);
        }
        printf("\nFINAL BIMODAL CONTENTS");
        for(int i=0;i<counter1->total_rows;i++)
        {
            printf("\n%d  %d",i,counter1->counter[i]);
        }
    }

    
    return 0;
}
