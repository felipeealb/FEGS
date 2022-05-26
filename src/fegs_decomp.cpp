//#include <iostream>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <time.h>
#include <sys/time.h>
#include <string>
#include <filesystem>
#include <list>
#include<queue>

#include<fstream>
using namespace std;

// ********** VARIAVEIS GLOBAIS **********
int num_vertices;/* numero de vertices/individuos do grafo*/
int num_skills;/*numero de habilidades*/
int num_teams; /* quantidade de total de equipes/projetos*/

int **G; /*matriz |n|x|n| grafo*/
int **Graph_Zuv; /*matriz |n|x|n| grafo ponderado pelos caminhos minimos positivos*/
int **K; /*matriz |n|x|f| de pessoa - habilidade*/
int **R; /* matriz demanda |m|x|f|: quantidade de cada pessoa com habilidade k_a em cada equipe*/

char instanceG[50];
char instanceKR[50];

char CURRENT_DIR[500];

//double runTime;
#define BILLION 1000000000L


//typedef
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
typedef IloArray<BoolVar3Matrix> BoolVar4Matrix;
//--------------------------------------------------------------------------------------------------------------
double get_wall_time(){
    struct timeval time;
    if(gettimeofday(&time,nullptr)){
        // HANDLE ERROR
        return 0;
    }else{
        return static_cast<double>(time.tv_sec) + static_cast<double>(time.tv_usec*0.000001); //microsegundos
    }
}

//--------------------------------------------read data------------------------------------------------------------------
//load matrix
void load_Graph(char arq[], bool show, string instance_type){

    FILE *inst = fopen(arq,"r");
    if(!inst)
        printf("erro na leitura do arquivo!");

    int edge_weight;
    fscanf(inst,"%d",&num_vertices);
    G = (int**)(malloc(num_vertices*sizeof(int*)));

    for (int u=0;u<num_vertices;u++)
    {
        G[u] = (int*)malloc(num_vertices*sizeof(int));

        for (int v=0;v<num_vertices;v++)
        {
            fscanf(inst,"%d",&edge_weight);
            // G[u][v] = edge_weight;
            if(u == v){
                G[u][v] = 1;
            }else{
                if(edge_weight>0)
                    G[u][v] = 1;
                else if(edge_weight<0)
                    G[u][v] = -1;
                else if(edge_weight==0)
                    G[u][v] = 0;
            }
        }

    }
    fclose(inst);

    if (instance_type == "random"){
        // cout << "random" << endl;
        for (int u = 0; u < num_vertices; u++){
            for (int v = u+1; v < num_vertices; v++){
                if (G[u][v] == 1 && G[v][u] != 1){
                    G[v][u] = 1;
                }else if (G[v][u] == 1 && G[u][v] != 1){
                    G[u][v] = 1;
                }else if (G[u][v] == -1 && G[v][u] != -1){
                    G[v][u] = -1;
                }else if (G[v][u] == -1 && G[u][v] != -1){
                    G[u][v] = -1;
                }
            }
        }
    }


    if(show){
        cout << endl << "Numero de individuos: " << num_vertices << endl;
        cout << " Adjacency Matrix \n";
        for (int u = 0; u < num_vertices; u++){
            for (int v = 0; v < num_vertices; v++)
                cout <<  G[u][v] << " ";
            cout<< "\n";
        }
    }
}

//load individuals skills
void load_K(char arq[], bool show){

    FILE *inst = fopen(arq,"r");
    if(!inst)
        printf("erro na leitura do arquivo!");

    int value;
    fscanf(inst,"%d",&num_skills);
    K = (int**)(malloc(num_vertices*sizeof(int*)));
    for (int u=0;u<num_vertices;u++){
        K[u] = (int*)malloc(num_skills*sizeof(int));

        for (int s=0;s<num_skills;s++){
            fscanf(inst,"%d",&value);
            K[u][s] = value;
        }
    }
    if(show){
        printf("\nHabilidades = %d \n", num_skills);
        for (int u=0;u<num_vertices;u++){
            for (int s=0;s<num_skills;s++)
                printf("%d ", K[u][s]);
            printf("\n");
        }
    }

    fclose(inst);
}

//load the skills necessary in projects/teams
void load_R(char arq[], bool show){

    FILE *inst = fopen(arq,"r");
    if(!inst)
        printf("erro na leitura do arquivo!");

    int value;
    fscanf(inst,"%d",&num_teams);
    R = (int**)(malloc(num_teams*sizeof(int*)));

    for (int j=0;j<num_teams;j++){
        R[j] = (int*)malloc(num_skills*sizeof(int));

        for (int s=0;s<num_skills;s++){
            fscanf(inst,"%d",&value);
            R[j][s] = value;
        }
    }
    if(show){
        printf("\nEquipes = %d\n", num_teams);

        for (int j=0;j<num_teams;j++){
            for (int s=0;s<num_skills;s++)
                printf("%d ", R[j][s]);
            printf("\n");
        }
    }

    fclose(inst);
}


static void
instanceReader(bool show, int vert, int graph_class, int class_type, string instance_type){
    char arq1[600]; char arq2[600]; char arq3[600];
    char arq_base[500];
    sprintf(arq_base, "%s/instances_COR", CURRENT_DIR);

    if (instance_type == "random"){
        sprintf(instanceG, "%dverticesS%d", vert, graph_class);
    }
    else if (instance_type == "bitcoinotc" || instance_type == "epinions"){
        sprintf(instanceG, "%dvertices_%s_S%d", vert, instance_type.c_str(), graph_class);

    }else{
        cout << "[INFO] not a valid instance type" << endl;
    }


    sprintf(arq1, "%s/%dVertices/%s.txt", arq_base, vert, instanceG);
    cout << "[INFO] Read graph: "<< instanceG << endl;

    sprintf(instanceKR, "%dVertices/class1/%d", vert, class_type);
    sprintf(arq2, "%s/%s/K.txt", arq_base, instanceKR);
    sprintf(arq3, "%s/%s/R.txt", arq_base, instanceKR);
    cout << "[INFO] Instance class type: "<< instanceKR << endl;


    load_Graph(arq1, show, instance_type); load_K(arq2, show); load_R(arq3, show);

}


static void
    instanceReader (bool show, int vert, int graph_class, int class_type, string instance_type);

static void
    allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y);

static void
    allocVars_uv (IloEnv env, BoolVarMatrix f, IloIntVar lambda);

static void
    createModel_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, int u, int v);

static void
   createModel_Decomp (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

static void
   objFunction_uv (IloModel model, BoolVarMatrix f, int u, int v);

static void
   objFunction_Decomp (IloModel model, BoolVar3Matrix y);

static void
   rest1 (IloModel model, BoolVar3Matrix x);

static void
   rest2 (IloModel model, BoolVar3Matrix x);

static void
   rest3 (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

static void
   rest4_uv (IloModel model, BoolVarMatrix f, int u, int v);


static void
   rest5_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, int u, int v);


static void
    displaySolution(IloCplex& cplex,
                   BoolVar3Matrix x,
                   BoolVar3Matrix y);

static void
    displaySolution_uv(IloCplex& cplex,
                   BoolVar4Matrix f,
                   IntVarMatrix lambda);


static void
    saveSolution(IloCplex& cplex,
                   BoolVar3Matrix X,
                   BoolVar3Matrix Y,
                   BoolVar4Matrix F,
                   int class_type);
static void
    saveResults(IloCplex& cplex,
                    double timeF);


static void
    savePath(IloCplex& cplex,
                   BoolVar4Matrix F,
                   int class_type);

static void
    createGraph_Zuv();


static void
    testeGrafos();

void traverse(int u, bool visited[]) {
   visited[u] = true; //mark v as visited
   for(int v = 0; v<num_vertices; v++) {
      if(G[u][v]) {
         if(!visited[v])
            traverse(v, visited);
      }
   }
}
bool isConnected() {
   bool *vis = new bool[num_vertices];
   //for all vertex u as start point, check whether all nodes are visible or not
   for(int u = 0; u < num_vertices; u++) {
      for(int i = 0; i<num_vertices; i++)
         vis[i] = false; //initialize as no node is visited
         traverse(u, vis);
      for(int i = 0; i<num_vertices; i++) {
         if(!vis[i]) //if there is a node, not visited by traversal, graph is not connected
            return false;
      }
   }
   return true;
}

int main()
{
    getcwd(CURRENT_DIR, 500);
    char intances_types[3][12] =   {"random", "bitcoinotc", "epinions"};
    // graph 25 vertices epinions_S3: disconnected (9-18-22) 
    
    double cpu0_exec, cpu1_exec;
    cpu0_exec = get_wall_time();
    int vert = 25;
    for(int i=0; i<3; i++){
        for(int gclass = 1; gclass<4; gclass++){  // 1-3
            for(int ctype = 1; ctype < 7; ctype++){  //1-6

                IloEnv env;

                try {

                    instanceReader(0,vert, gclass,ctype, intances_types[i]);

                    double timeZ_start = get_wall_time();
                    createGraph_Zuv();
                    double timeZ_final = get_wall_time();
                    double timeGraph = timeZ_final - timeZ_start;
                    cout << "[INFO]: Create graph Zuv" << endl;
                    
                    // create ILP problem
                    IloModel model(env);
                    cout << "[INFO]: Create variables decomposed" << endl;
                    BoolVar3Matrix x(env,num_vertices);
                    BoolVar3Matrix y(env,num_vertices);
                    allocVars(env, x, y);
                    cout << "[INFO]: Create model decomposed" << endl;
                    createModel_Decomp(model, x, y);
                    IloCplex cplex;
                    cplex = IloCplex(model);
                    // cplex.exportModel("FEGS_Decomposed.lp");
                       
                    double cpu0, cpu1;
                    cplex.setParam(IloCplex::TiLim, 7200); // time limit 2h
                    cplex.setParam(IloCplex::TreLim, 7000); // memory limit 7GB
                    // cplex.setParam(IloCplex::WorkMem, 2000);

                    cpu0 = get_wall_time();
                    if (!cplex.solve()){
                        env.error() << "[INFO]: Failed to optimize ILP." << endl;
                        cout << "Solution status = " << cplex.getStatus()   << endl;
                        throw(-1);
                    }
                    cpu1 = get_wall_time();
                    double runTime = (cpu1 - cpu0) + timeGraph;
                    

                    cout << "Solution status = " << cplex.getStatus()   << endl;
                    cout << "Solution value  = " << cplex.getObjValue() << endl;
                    cout << "cplex time: " << cplex.getTime() << endl;
                    cout << "run time: " << runTime << endl;

                    // displaySolution(cplex,x,y);
                    // saveSolution(cplex,x,y,f,ctype);
                    saveResults(cplex,runTime);


                }catch (IloException& e) {
                    cerr << "Concert exception caught: " << e << endl;
                }
                catch (...) {
                    cerr << "Unknown exception caught" << endl;
                }

            env.end();

            }

        }
    }

    cout << "--------------------------------------------------------" << endl;
    cpu1_exec = get_wall_time();
    cout << "Total time: " << cpu1_exec - cpu0_exec << endl;

    return 0;
}


static void
allocVars_uv (IloEnv env, BoolVarMatrix f, IloIntVar lambda){

    // var f(u,v)[p][q]: if arc pq is used in the flow associate to path(u,v)
    for (int p=0; p<num_vertices; p++)
        f[p] = IloBoolVarArray(env, num_vertices);
    
    

}

static void
allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y){

    // var x[u][j][s] : if worker u is in projecet j with skill s
    for (int u=0; u<num_vertices; u++){
        x[u] = BoolVarMatrix(env, num_teams);
        for (int j = 0; j < num_teams; j++)
            x[u][j] = IloBoolVarArray(env, num_skills);
    }

    // var y[u][v][j] : if worker u and v is in the same project j
    for (int u=0; u<num_vertices; u++){
        y[u] = BoolVarMatrix(env, num_vertices);
        for (int v = u+1; v < num_vertices; v++) // u<v
            y[u][v] = IloBoolVarArray(env, num_teams);
    }


}


static void
createModel_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, int u, int v){

    // add var f(u,v)[p][q]
    if(G[u][v] < 1){ // and u,v not in E{+}
        for(int p=0; p<num_vertices; p++)
            for(int q=0; q<num_vertices; q++) // for all arc (p,q)
                if(p!=q && G[p][q] != 0){ // must have an edge(p,q) to have an arc (p,q)
                    char name[50];
                    sprintf(name, "F(%d,%d)%d%d", u+1, v+1,p+1, q+1);
                    f[p][q].setName(name);
                    model.add(f[p][q]);
                }
    }

    // add var lambda(u,v)
    if(G[u][v] < 1){ // and u,v not in E{+}
        char name[20];
        sprintf(name, "lambd(%d,%d)", u+1, v+1);
        lambda.setName(name);
        model.add(lambda);
    }

    objFunction_uv (model,f,u,v);
    rest4_uv(model,f, u, v); // flow const 
    rest5_uv(model,f,lambda, u,v); // every path with neg edges is pair



}


static void
createGraph_Zuv(){

    Graph_Zuv = (int**)(malloc(num_vertices*sizeof(int*)));
    for (int u=0;u<num_vertices;u++)
        Graph_Zuv[u] = (int*)malloc(num_vertices*sizeof(int));

    for(int u=0; u<num_vertices;u++)
        for(int v=0; v<num_vertices; v++)
            if(u==v)
                Graph_Zuv[u][v] = 1;

    for(int u=0; u<num_vertices;u++)
        for(int v=u+1; v<num_vertices; v++)
            if(G[u][v] < 1 ){

                IloEnv env;
                try{
                    IloModel model(env);
                    // cout << "[INFO]: Create variables" << endl;
                    BoolVarMatrix f(env,num_vertices);
                    IloIntVar lambda(env);
                    allocVars_uv(env,f,lambda);
                    createModel_uv(model,f,lambda,u,v);        // can optmize var f, put all dif u,v = zero or better use 2 indices
                    IloCplex cplex;
                    cplex = IloCplex(model);
                    cplex.setOut(env.getNullStream());
                    if (!cplex.solve()){
                        // env.error() << "[INFO]: Failed to optimize ILP." << endl;
                        // cout << "Solution status = " << cplex.getStatus()   << endl;
                        Graph_Zuv[u][v] = 1000000;
                        Graph_Zuv[v][u] = 1000000;
                        throw(-1);
                    }

                    // cout << "Solution status = " << cplex.getStatus()   << endl;
                    // cout << "Solution value  = " << cplex.getObjValue() << endl;

                    Graph_Zuv[u][v] = int(cplex.getObjValue());
                    Graph_Zuv[v][u] = int(cplex.getObjValue());

                    // cout << "-------------------------------------------------- \n\n" << endl;
                }catch (IloException& e) {
                        cerr << "Concert exception caught: " << e << endl;
                    }
                    catch (...) {
                        // cerr << "Unknown exception caught" << endl;
                    }

                env.end();


            }else{
                Graph_Zuv[u][v] = 1;
                Graph_Zuv[v][u] = 1;
            }
    // cout << "-------------------------------" << endl;
    // cout << " Adjacency Matrix Zuv\n";
    // for (int u = 0; u < num_vertices; u++){
    //     for (int v = 0; v < num_vertices; v++)
    //         cout <<  Graph_Zuv[u][v] << " ";
    //     cout<< "\n";
    // }

}


static void
createModel_Decomp (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y){

    // add var x[u][j][s] in model
    for(int j=0; j<num_teams; j++) // for all project j 
        for(int u=0; u< num_vertices; u++) // for all worker u
            for(int s=0; s<num_skills; s++){ // 
                if(K[u][s] > 0){ // for all skills: s(u)
                    char name[30];
                    sprintf(name, "x%d%d%d", u+1, j+1,s+1);
                    x[u][j][s].setName(name);
                    model.add(x[u][j][s]);
                }
            }

    // add var y[u][v][j]
    for(int u=0; u<num_vertices; u++) // for all u 
        for(int v = u+1; v<num_vertices; v++) // for all v: u<v
            for(int j = 0; j<num_teams; j++){ // for all project j
                char name[30];
                sprintf(name, "y%d%d%d", u+1, v+1,j+1);
                y[u][v][j].setName(name);
                model.add(y[u][v][j]);
            }


    objFunction_Decomp(model,y);
    rest1(model, x); // max one team with one skill
    rest2(model, x); // min skill s per team j
    rest3(model,x,y); // linearization y with x

}


static void
objFunction_uv (IloModel model, BoolVarMatrix f, int u, int v){
    IloEnv env = model.getEnv();
    IloExpr objExpr(env);
    for(int p=0; p<num_vertices; p++)
        for(int q=0; q<num_vertices; q++)
            if(G[u][v] < 1 && p!= q && G[p][q] != 0){ // (u,v) not in E+ and (p,q) is a possible arc between (u,v)-path
                objExpr += f[p][q];
            }

    IloObjective obj = IloMinimize(env, objExpr);
    model.add(obj);
    objExpr.end();

}

static void
objFunction_Decomp (IloModel model, BoolVar3Matrix y){
    IloEnv env = model.getEnv();

    IloExpr objExpr(env);

    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++)
            for(int j = 0; j<num_teams; j++){
                objExpr += Graph_Zuv[u][v]*y[u][v][j];  // (u,v) in E+
            }

    IloObjective obj = IloMinimize(env, objExpr);
    model.add(obj);
    objExpr.end();
}

static void
rest1 (IloModel model, BoolVar3Matrix x){

    IloEnv env = model.getEnv();

    // each worker uses at most one skill 
    for (int u=0; u<num_vertices; u++) {
        IloExpr expr(env);
        for (int j = 0; j <num_teams; j++) {
            for(int s = 0; s<num_skills; s++){
                if(R[j][s]>0 && K[u][s]>0) // team j need skill s and indidual u have skill s
                    expr += x[u][j][s];
            }
        }
        model.add(expr <= 1);   // posso adcionar variavael boleana para saber sse sumConst está empty ou nao
        expr.end();
    }

}

static void
rest2 (IloModel model, BoolVar3Matrix x){
    IloEnv env = model.getEnv();

    // minimum number of individuals in project j with skill s
    for (int j=0; j <num_teams; j++) {
        for (int s=0; s <num_skills; s++) {
            IloExpr expr(env);
            if (R[j][s] > 0) { //for all skill s in project j
                for (int u=0; u < num_vertices; u++) {
                    if(K[u][s] > 0) // if worker u have skill s
                        expr += x[u][j][s];
                }
                model.add(expr >= R[j][s]);
                expr.end();
            }
        }
    }

}

static void
rest3 (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y){
    IloEnv env = model.getEnv();

    // set var y with x (linearization)
    for (int u=0; u<num_vertices; u++)
        for (int v = u+1; v < num_vertices; v++) // for all u,v: u<v
            for (int j = 0; j <num_teams; j++){
                IloExpr exprU(env);
                IloExpr exprV(env);
                for (int s = 0; s < num_skills; s++) {
                    if(K[u][s] > 0 && R[j][s] > 0) // individual u have skill s and team j need skill s
                        exprU += x[u][j][s];
                    if(K[v][s] > 0 && R[j][s] > 0)
                        exprV += x[v][j][s];
                }
                model.add(exprU + exprV - y[u][v][j] <= 1);
                model.add(y[u][v][j]  <= exprU);
                model.add(y[u][v][j] <= exprV);
                exprU.end(); exprV.end();
            }

}

static void
rest4_uv (IloModel model, BoolVarMatrix f, int u, int v){
    IloEnv env = model.getEnv();

    if(G[u][v] < 1){  // for all u<v and (u,v) not in E+
        for(int q = 0; q<num_vertices; q++){ // for fixed vertice q
            IloExpr sumArcs(env);
            for(int p = 0; p<num_vertices; p++) // for every arc (p,q) <=> (p,q) or (q,p) in E 
                if(p!=q && G[p][q] != 0){ //  (p,q) in E <=> (q,p) in E
                    sumArcs += f[p][q];
                    sumArcs += -f[q][p];
                }
            if(q == u){
                IloExpr expr(env);
                model.add(sumArcs == -1); expr.end();
            }else if(q == v){
                IloExpr expr(env);
                model.add(sumArcs == 1); expr.end();
            }else{
                model.add(sumArcs == 0);
            }
        }
    }

}

static void
rest5_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, int u, int v){
    IloEnv env = model.getEnv();
    
    if(G[u][v] < 1 ){ // (u,v) not in E+
        IloExpr expr(env);
        for(int p = 0; p<num_vertices; p++)
            for(int q = 0; q<num_vertices; q++)
                if(p!=q && G[p][q] < 0){ // sum of negatives arcs(A-)
                    expr += f[p][q];
                }
        model.add(expr - 2*lambda == 0); // pair number of negatives edges
        expr.end();
    }
        
}


static void
displaySolution(IloCplex& cplex,
               BoolVar3Matrix x,
               BoolVar3Matrix y){



    cout << "Solution status = " << cplex.getStatus()   << endl;
    cout << "Solution value  = " << cplex.getObjValue() << endl;
    cout << "Valores das Variaveis " << endl;

    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------" << endl;

    for(int j=0; j<num_teams; j++)
        for(int u=0; u<num_vertices; u++)
            for(int s=0; s<num_skills; s++){
                if(K[u][s] > 0){
                    if(cplex.getValue(x[u][j][s]) > 0){
                        cout << "\t " << "x"<< "["<< u+1  << "," << j+1 << "," << s+1 << "]" << " -> " << " 1 \n";
                    }
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++)
            for(int j=0; j<num_teams; j++){
                if(cplex.getValue(y[u][v][j]) > 0){
                    cout << "\t " << "y"<< "["<< u+1  << "," << v+1 << "," << j+1 << "]" << " -> " << " 1 \n";
                }
            }

}

static void
    displaySolution_uv(IloCplex& cplex,
                   BoolVar4Matrix f,
                   IntVarMatrix lambda){

    cout << "--------------------------------------------------------" << endl;

    for(int u = 0; u<num_vertices; u++)
        for(int v = 0; v<num_vertices; v++)
            if(u<v && G[u][v] < 1){
                for(int p = 0; p<num_vertices; p++)
                    for(int q = 0; q<num_vertices; q++)
                        if(p!= q && G[p][q] != 0){
                            if(cplex.getValue(f[u][v][p][q]) > 0){
                                cout << "\t " << "f"<< "["<< u+1  << "," << v+1 << "," << p+1 << ", " << q+1 << "]" << " -> " << " 1 \n";
                            }
                        }
            }


}

static void
    saveSolution(IloCplex& cplex,
                   BoolVar3Matrix x,
                   BoolVar3Matrix y,
                   BoolVar4Matrix f,
                   int class_type){


    char arq[1000];
    char arqv_instance[50];
    sprintf(arqv_instance, "%s_%d", instanceG,class_type);
    sprintf(arq, "%s/results/%s.txt",CURRENT_DIR, arqv_instance);

    FILE *file = fopen(arq, "w");
    if (file == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(file, "Solution value = %.1f \n", cplex.getObjValue());

    for(int j=0; j<num_teams; j++)
        for(int u=0; u<num_vertices; u++)
            for(int s=0; s<num_skills; s++){
                if(K[u][s] > 0){
                    if(cplex.getValue(x[u][j][s]) > 0){
                        fprintf(file, "x[%d,%d,%d] = 1 \n", u+1, j+1, s+1);
                    }
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++)
            for(int j=0; j<num_teams; j++){
                if(cplex.getValue(y[u][v][j]) > 0){
                    fprintf(file, "y[%d,%d,%d] = 1 \n", u+1, v+1, j+1);
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u = 0; u<num_vertices; u++)
        for(int v = 0; v<num_vertices; v++)
            if(u<v && G[u][v] < 1){
                for(int p = 0; p<num_vertices; p++)
                    for(int q = 0; q<num_vertices; q++)
                        if(p!= q && G[p][q] != 0){
                            if(cplex.getValue(f[u][v][p][q]) > 0){
                                fprintf(file, "f[%d,%d,%d,%d] = 1 \n", u+1, v+1, p+1,q+1);
                            }
                        }
            }


    fclose(file);

}


static void
    saveResults(IloCplex& cplex,
                   double timeF){

    char arq[600];
    sprintf(arq, "%s/results/%d_Vertices_17-05-2022_Decomp.ods",CURRENT_DIR, num_vertices);
    
    ofstream outputTable;
    outputTable.open(arq,ios:: app);
    if(outputTable.is_open()){

        outputTable << instanceG << ";"; // grafo instancia
        outputTable << instanceKR << ";"; // KR instancia
        outputTable << num_vertices << ";";   // numero de vertices
        outputTable << num_skills << ";";   // qtd de hab
        outputTable << num_teams << ";";   // qtd de proj
        outputTable << cplex.getStatus() << ";"; // Status cplex
        outputTable << cplex.getObjValue() << ";"; // valor fo
        outputTable << cplex.getNnodes() << ";"; // num nos
        outputTable << cplex.getMIPRelativeGap() <<";"; // gap
        outputTable << timeF <<  ";"; // tempo execucao flp
        outputTable << cplex.getTime() <<  ";"; // tempo execucao cplex
        outputTable << " \n ";


    }

    outputTable.close();

}


static void
    testeGrafos(){

    char intances_types[3][12] =   {"random", "bitcoinotc", "epinions"};
    
    cout << "-------------------------------------------------" << endl;
    int vert = 25;
    for(int i=0; i<3; i++){
        for(int gclass = 1; gclass<4; gclass++){  // 1-3
            
            instanceReader(0,vert, gclass,1, intances_types[i]);
            
            if(isConnected())
                cout << "The Graph is connected." << endl;
            else
                cout << "The Graph is not connected." << endl;


            cout << "-------------------------------------------" << endl;
            
        }
    }
    
}

