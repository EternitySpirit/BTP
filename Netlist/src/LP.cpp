#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <stack>
// #include "ortools/linear_solver/linear_solver.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"

#define MAXFANIN 5
#define MAXLEVEL 500
#define MAXPI 100
#define MAXGATES 10000

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
std::map<int,std::vector<int>> adj_list;

void parse_gates();
void compute_asap_level();
void compute_alap_level();
void compute_list_level();
void print_gates();
void print_stat();
int count_v(const char* line);
void get_vars(const char* line, int varid[], int nv);
bool all_alap_labeled();
void update_alap(int gate_index, int output);

using namespace operations_research;
using namespace sat;

void MinimizeFinalTimeStepWithConstraints(int numGates, int M, int maxTimeSteps, const std::vector<table_gate>& gates) {
    CpModelBuilder model;

    // Define variables: x[i][t] is True if gate i is assigned to time step t.
    std::vector<std::vector<BoolVar>> x(numGates, std::vector<BoolVar>(maxTimeSteps));
    for (int i = 0; i < numGates; ++i) {
        for (int t = 0; t < maxTimeSteps; ++t) {
            x[i][t] = model.NewBoolVar();
        }
    }

    // Define the final time variable: an integer variable representing the final time step.
    const Domain time_domain(0, maxTimeSteps - 1);
    const IntVar final_time = model.NewIntVar(time_domain);

    // Objective: Minimize final_time.
    model.Minimize(final_time);

    // Each gate should be scheduled exactly once.
    for (int i = 0; i < numGates; ++i) {
        LinearExpr sum;
        for (int t = 0; t < maxTimeSteps; ++t) {
            sum += x[i][t];
        }
        model.AddEquality(sum, 1);
    }

    // Limit the number of gates that can be scheduled at each time step to M.
    for (int t = 0; t < maxTimeSteps; ++t) {
        LinearExpr row_sum;
        for (int i = 0; i < numGates; ++i) {
            row_sum += x[i][t];
        }
        model.AddLessOrEqual(row_sum, M);
    }

        // Dependency constraints: if gate i depends on predecessor gate, it should be scheduled afterward.
       for (int i = 0; i < numGates; ++i) {
        for (int j = 0; j < gates[i].fanin; ++j) {
            int predecessor = gates[i].input[j];
            if (predecessor >= 0) {
                for (int t = 0; t < maxTimeSteps - 1; ++t) {
                  //Constraint to ensure i starts only after predecessor is completed
                  for (int s = 0; s<= t; s++)
                  {
                    model.AddImplication(x[predecessor][t], x[i][s].Not());
                  }
                  
                }
            }
        }
    }
   // Find leaves using DFS and apply the final_time constraint only to them.
    std::vector<bool> visited(numGates, false);
    std::vector<int> leaves;

    for (int i = 0; i < numGates; ++i) {
        if (!visited[i]) {
            std::stack<int> stack;
            stack.push(i);
            visited[i] = true;

            while (!stack.empty()) {
                int node = stack.top();
                stack.pop();

                // If the gate has no outputs, it's a leaf.
                if (gates[node].output.empty()) {
                    leaves.push_back(node);
                }

                // Visit all children of the current node.
                for (int child : gates[node].output) {
                    if (!visited[child]) {
                        stack.push(child);
                        visited[child] = true;
                    }
                }
            }
        }
    }

    // Apply final_time constraint only to leaves.
    for (int leaf : leaves) {
        for (int t = 0; t < maxTimeSteps; ++t) {
            model.AddGreaterOrEqual(final_time, t).OnlyEnforceIf(x[leaf][t]);
        }
    }

   // Solving the model.
    Model cp_model;
    CpSolverResponse response = SolveCpModel(model.Build(), &cp_model);

    if (response.status() == CpSolverStatus::OPTIMAL) {
        std::cout << "Optimal final time: " << SolutionIntegerValue(response, final_time) << std::endl;
        for (int i = 0; i < numGates; ++i) {
            for (int t = 0; t < maxTimeSteps; ++t) {
                if (SolutionBooleanValue(response, x[i][t])) {
                    std::cout << "Gate " << i << " assigned to time step " << t << std::endl;
                }
            }
        }
    } else {
        std::cerr << "No optimal solution found!" << std::endl;
    }
}

// Function to find leaf nodes using DFS
void find_leaves_with_final_time_constraint(int *leaves, int *leaf_count) {
    bool visited[MAXGATES] = {false};  // Track visited nodes
    int stack[MAXGATES];
    int stack_top = -1;
    *leaf_count = 0;

    for (int i = 0; i < ng; i++) {
        if (!visited[i]) {
            // Initialize DFS stack with the current node
            stack[++stack_top] = i;
            visited[i] = true;

            while (stack_top >= 0) {
                int node = stack[stack_top--];  // Pop from stack

                // If the gate has no outputs, it's a leaf node
                if (gates[node].output_count == 0) {
                    leaves[(*leaf_count)++] = node;
                }

                // Push all output gates (children) of the current node to stack
                for (int j = 0; j < gates[node].output_count; j++) {
                    int child = gates[node].output_gates[j];
                    if (!visited[child]) {
                        stack[++stack_top] = child;
                        visited[child] = true;
                    }
                }
            }
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
  for x20, the integer 5020 is stored in varid[].
*******************************************************/
void get_vars(char x[], int varid[], int nv)
{
  int i = 0, tval, count = 0, nx_flag;

  do {
    // Skip to the next digit in x[]
    while (!isdigit(x[i])) i++;

    // Check if the variable starts with 'x'
    if (x[i-1] == 'x') nx_flag = 1;
    else nx_flag = 0;

    // Convert number to integer
    tval = 0;
    while (isdigit(x[i])) {
      tval = tval * 10 + (x[i] - '0');
      i++;
    }

    // If it starts with 'x', add 5000
    if (nx_flag == 1) tval = 5000 + tval;
    varid[count++] = tval;

  } while (count < nv);
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
    for (int j = 0; j < gates[ng-1].fanin; j++)
    {
      if (varid[j]<5000)
      {
        adj_list[ng].push_back(varid[j])
      }

    }
  }
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
bool all_alap_labeled()
{
  for (int i = 0; i < ng; i++) {
    if (gates[i].alap_level == -1)  // If any gate is not labeled
      return false;
  }
  return true;  // All gates are labeled
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


int main(int argc, char *argv[]) {
    const char *file_name = "./src/net1.txt";
    netlist = fopen(file_name, "r");

    if (!netlist) {
        std::cerr << "Error opening netlist file." << std::endl;
        return 1;
    }

    FILE *output_file = freopen("f1.txt", "w", stdout);

    parse_gates();
    compute_asap_level();
    compute_alap_level();

    int M = 1;
    int maxTimeSteps = 4;

    std::cout << "\n************ BENCHMARK: " << argv[1] << " *************\n";
    //print_gates();
    //print_stat();

    operations_research::MinimizeFinalTimeStepWithConstraints(ng, M, maxTimeSteps);

    if (output_file) {
        fclose(output_file);
    }
    fclose(netlist);
    return 0;
}

