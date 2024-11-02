#include <iostream>
#include <vector>
#include <limits>
#include "ortools/linear_solver/linear_solver.h"

#define MAXGATES 5000
#define MAXFANIN 5
#define MAXLEVEL 500
#define MAXPI 100

struct table_gate {
    int fanin;
    int input[MAXFANIN];
    int out;
    int asap_level;
    int alap_level;
    int mobility;
    int list_level;
};

// Global Variables
table_gate gates[MAXGATES];
int primary_inputs[MAXPI]; 
int nPI = 0;
int ng = 0;
int tempvar = 1;
FILE *netlist;

// Forward declarations for existing C functions
void parse_gates();
void compute_asap_level();
void compute_alap_level();
void compute_list_level();
void print_gates();
void print_stat();

namespace operations_research {

void MinimizeFinalTimeStepWithConstraints(int numGates, int M, int maxTimeSteps) {
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("CBC"));
    if (!solver) {
        std::cerr << "Solver not available." << std::endl;
        return;
    }

    // Variables: x[i][t] = 1 if gate i is scheduled at time t, else 0
    std::vector<std::vector<MPVariable*>> x(numGates, std::vector<MPVariable*>(maxTimeSteps));
    for (int i = 0; i < numGates; ++i) {
        for (int t = 0; t < maxTimeSteps; ++t) {
            x[i][t] = solver->MakeBoolVar("x_" + std::to_string(i) + "_" + std::to_string(t));
        }
    }

    // Final time variable for objective minimization
    MPVariable* final_time = solver->MakeIntVar(0.0, maxTimeSteps - 1, "FinalTime");

    // Objective: Minimize the final time step
    solver->Minimize(final_time);

    // 1. Task Assignment Constraints
    // Ensure that for each time step, the number of gates assigned does not exceed M
    for (int t = 0; t < maxTimeSteps; ++t) {
        MPConstraint* row_constraint = solver->MakeRowConstraint(0, M, "RowLimit_" + std::to_string(t));
        for (int i = 0; i < numGates; ++i) {
            row_constraint->SetCoefficient(x[i][t], 1);
        }
    }

    // 2. Dependency Constraints
    // Ensure each gate starts only after all its predecessors are completed
    for (int i = 0; i < numGates; ++i) {
        for (int j = 0; j < gates[i].fanin; ++j) {
            int predecessor = gates[i].input[j];
            if (predecessor >= 0) {  // Only apply if dependency is valid
                for (int t = 0; t < maxTimeSteps - 1; ++t) {
                    // Constraint to ensure i starts only after predecessor is completed
                    for (int s = 0; s <= t; s++)
                    {
                        x[predecessor][t]>x[i][s]
                    }
                    
                }
            }
        }
    }

    // Solve the ILP
    const MPSolver::ResultStatus result_status = solver->Solve();
    if (result_status == MPSolver::OPTIMAL) {
        std::cout << "Optimal final time: " << final_time->solution_value() << std::endl;
        for (int i = 0; i < numGates; ++i) {
            for (int t = 0; t < maxTimeSteps; ++t) {
                if (x[i][t]->solution_value() == 1) {
                    std::cout << "Gate " << i << " assigned to time step " << t << std::endl;
                }
            }
        }
    } else {
        std::cerr << "No optimal solution found!" << std::endl;
    }
}

}  // namespace operations_research

int main(int argc, char *argv[]) {
    const char *file_name = "BENCH/Bench/xor5_d.txt";
    netlist = fopen(file_name, "r");

    if (!netlist) {
        std::cerr << "Error opening netlist file." << std::endl;
        return 1;
    }

    FILE *output_file = freopen("f1.txt", "w", stdout);

    parse_gates();
    compute_asap_level();
    compute_alap_level();

    int M = 4;           // Crossbar array max rows per level
    int maxTimeSteps = 100; // Adjust maximum time steps as per circuit complexity

    std::cout << "\n************ BENCHMARK: " << argv[1] << " *************\n";
    print_gates();
    print_stat();

    operations_research::MinimizeFinalTimeStepWithConstraints(ng, M, maxTimeSteps);

    fclose(netlist);
    fclose(output_file);

    return 0;
}

void parse_gates()
{
  char line[50];
  int  varid[10], nvar;


  for (;;)
  {
    fgets (line, sizeof(line), netlist);

    if (line[0] == '.') break;

    nvar = count_v (line);   // count how many variables (i.e. n or x)
    get_vars (line, varid, nvar);

    switch (nvar)
    {
      case 2: // NOT gate
              gates[ng].fanin = 1;
  	      gates[ng].out = varid[0];
	      gates[ng].input[0] = varid[1];
	      ng++;
	      break;

      case 3: // 2-input NOR
              gates[ng].fanin = 2;
	      gates[ng].out = varid[0];
	      gates[ng].input[0] = varid[1];
	      gates[ng].input[1] = varid[2];
	      ng++;
	      break;

      case 4: // Two 2-input NOR gates in cascade
              gates[ng].fanin = 2;
	      gates[ng].out = -(tempvar++);
	      gates[ng].input[0] = varid[2];
	      gates[ng].input[1] = varid[3];
	      ng++;
              gates[ng].fanin = 2;
	      gates[ng].out = varid[0];
	      gates[ng].input[0] = varid[1];
	      gates[ng].input[1] = -(tempvar-1);
	      ng++;
	      break;

      case 5: // Three 2-input NOR gates in two levels
              gates[ng].fanin = 2;
	      gates[ng].out = -(tempvar++);
	      gates[ng].input[0] = varid[1];
	      gates[ng].input[1] = varid[2];
	      ng++;
              gates[ng].fanin = 2;
	      gates[ng].out = -(tempvar++);
	      gates[ng].input[0] = varid[3];
	      gates[ng].input[1] = varid[4];
	      ng++;
              gates[ng].fanin = 2;
	      gates[ng].out = varid[0];
	      gates[ng].input[0] = -(tempvar-2);
	      gates[ng].input[1] = -(tempvar-1);
	      ng++;
	      break;

      default: printf ("** Invalid variables on line: %s\n", line);
               exit;
    }
  }
}



/******************************************************
    Count the number of variables in line in x[] 
******************************************************/
int count_v (char x[])
{
  int i, k, count;

  k = strlen(x);
  count = 0;

  for (i=0; i<k; i++)
    if ((x[i]=='n') || (x[i]=='x')) count++;

  return count;
}


/*******************************************************
  Read the string x[], and store the ids of all the
  variables in the array varid[]. If a variable starts
  with 'x', 5000 is added to its id.
  For n20, the integer 20 is stored in varid[];
  for x20, the integer 5020 is stored in varid[]
*******************************************************/
get_vars (char x[], int varid[], int nv)
{
  int i, tval, count, nx_flag;

  i=0; count=0;
  do
    {
      while (!isdigit(x[i])) i++;
      
      if (x[i-1] == 'x') nx_flag=1;  else nx_flag=0;

      tval = 0;
      while (isdigit(x[i])) {
        tval = tval*10 + (x[i] - '0');
	i++;
      } 

      if (nx_flag == 1) tval = 5000 + tval;
      varid[count++] = tval;
    }
  while (count < nv); 
}


/*************************************************
  Tests whether a givel line is a PI. If so, it
  returns 1; else returns 0.
*************************************************/
int isPI (int lineid)
{
  int i;

  for (i=0; i<ng; i++)
    if (gates[i].out == lineid)  // Not a PI
      return 0;        

  return 1;
}


/*************************************************
  Tests whether a givel line is a PO. If so, it
  returns 1; else returns 0.
*************************************************/
int isPO (int lineid)
{
  int i, j;

  for (i=0; i<ng; i++)
  {
    for (j=0; j<gates[i].fanin; j++)
      if (gates[i].input[j] == lineid)  // Not a PO
        return 0;        
  }

  return 1;
}


/**************************************************
  Return the level of a given gate/PI identified
  by an id (integer); returns -1 if the id has 
  not yet been assigned a level.
**************************************************/
int get_asap_level (int lineid)
{
  int i;

  if (isPI(lineid))  return 0;

  for (i=0; i<ng; i++)
    if (gates[i].out == lineid)
      return (gates[i].asap_level);

  printf ("\n*** ERROR in get_asap_level ***\n");
  exit;
}



/**************************************************
  Return the level of a given gate/PI identified
  by an id (integer); returns -1 if the id has 
  not yet been assigned a level.
**************************************************/
int get_list_level (int lineid)
{
  int i;

  if (isPI(lineid))  return 0;

  for (i=0; i<ng; i++)
    if (gates[i].out == lineid)
      return (gates[i].list_level);

  printf ("\n*** ERROR in get_list_level ***\n");
  exit;
}


/**************************************************
  Return the level of a given gate/PI identified
  by an id (integer); returns -1 if the id has 
  not yet been assigned a level.
**************************************************/
int get_alap_level (int lineid)
{
  int i;

  if (isPO(lineid))  return 0;

  for (i=0; i<ng; i++)
    if (gates[i].out == lineid)
      return (gates[i].alap_level);

  printf ("\n*** ERROR in get_alap_level ***\n");
  exit;
}



/*****************************************************
  Fill up the gate levels for ASAP scheduling
*****************************************************/
void compute_asap_level()
{
  int i, t1, t2, max, count=0;

  while (count < ng)  // Till all gates are levelized
  {
    for (i=0; i<ng; i++)
    {
      if (gates[i].fanin == 1)  // NOT gate
      {
        t1 = get_asap_level(gates[i].input[0]);
        if (t1 != -1)
	{
	  gates[i].asap_level = t1 + 1;
	  count++;
	}
      }


      if (gates[i].fanin == 2)
      {
        t1 = get_asap_level(gates[i].input[0]);
        t2 = get_asap_level(gates[i].input[1]);

	if (t1>t2) max=t1;  else max=t2;

	if ((t1 != -1) && (t2 != -1))
	{
	  gates[i].asap_level = max + 1;
	  count++;
	}
      }
    }
  }
}


/*****************************************************
  Fill up the gate levels for ALAP scheduling
*****************************************************/
void compute_alap_level()
{
  int i, iter, maxl;

  for (i=0; i<ng; i++)
    if (isPO(gates[i].out))
       gates[i].alap_level = 1;

  while (! all_alap_labeled())  // Till all gates are levelized
  {
    for (i=0; i<ng; i++)
      update_alap (i, gates[i].out);
  }

  for (iter=0; iter<10; iter++)
    for (i=0; i<ng; i++)
      update_alap (i, gates[i].out);

  // Correct ALAP levels
  maxl = 0;
  for (i=0; i<ng; i++)
    if (gates[i].alap_level > maxl)
      maxl = gates[i].alap_level;
printf ("\n*** %d ***", maxl);
  for (i=0; i<ng; i++)
    gates[i].alap_level = maxl - gates[i].alap_level + 1;

}


/****************************************************
  Test whether all the gates have been labeled using
  ALAP scheduling.
****************************************************/
int all_alap_labeled()
{
  int i, flag;

  flag=0;

  for (i=0; i<ng; i++)
    if (gates[i].alap_level == -1)  return 0;

  return 1;
}


/***************************************************
  Update the labels of the gates during creation
  of the ALAP schedule.
***************************************************/
void update_alap (int index, int lineid)
{
  int j, k;

  for (j=0; j<ng; j++)
  {
    for (k=0; k<gates[j].fanin; k++)
    {
      if (lineid == gates[j].input[k])
        if (gates[j].alap_level != -1)
	  if (gates[index].alap_level <= gates[j].alap_level)
	    gates[index].alap_level = gates[j].alap_level + 1;
    }
  }
}

