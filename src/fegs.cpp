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
void load_Graph(char arq[], bool show){

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
            if(u == v){
                G[u][v] = 1;
            }else{
                if(edge_weight>0)
                    G[u][v] = 1;
                else if(edge_weight<0)
                    G[u][v] = -1;
            }
        }

    }
    fclose(inst);

//     Triangular Inferior
//    for (int u = 0; u < num_vertices; u++){
//        for (int v = u+1; v < num_vertices; v++)
//            G[u][v] = G[v][u];
//    }

//     Triangular Superior
    // for (int u = 0; u < num_vertices; u++){
    //     for (int v = u+1; v < num_vertices; v++)
    //         G[v][u] = G[u][v];
    // }

//     Dense
   for (int u = 0; u < num_vertices; u++){
       for (int v = u+1; v < num_vertices; v++){
           if (G[u][v] == 1 || G[v][u] == 0){
               G[v][u] = 1;
           }else if (G[v][u] == 1 || G[u][v] == 0){
               G[u][v] = 1;
           }else if (G[u][v] == -1 || G[v][u] == 0){
               G[v][u] = -1;
           }else if (G[v][u] == -1 || G[u][v] == 0){
               G[u][v] = -1;
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


    load_Graph(arq1, show); load_K(arq2, show); load_R(arq3, show);

}


static void
    instanceReader (bool show, int vert, int graph_class, int class_type, string instance_type);

static void
   allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda);

static void
   createModel (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda);

static void
   objFunction (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

static void
   rest1 (IloModel model, BoolVar3Matrix x);

static void
   rest2 (IloModel model, BoolVar3Matrix x);

static void
   rest3 (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

static void
   rest4 (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

static void
   rest5 (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

static void
   rest6 (IloModel model, BoolVar4Matrix f, IntVarMatrix lambda);


static void
    displaySolution(IloCplex& cplex,
                   BoolVar3Matrix x,
                   BoolVar3Matrix y,
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
    testeGrafos();

int main()
{
    getcwd(CURRENT_DIR, 500);
    char intances_types[3][12] =   { "random", "bitcoinotc", "epinions"};

    double cpu0_exec, cpu1_exec;
    // instanceReader(1,5, 1,1, "random");
    cpu0_exec = get_wall_time();
    int vert = 4;
    for(int i=0; i<1; i++){
        for(int gclass = 1; gclass<2; gclass++){
            for(int ctype = 1; ctype < 2; ctype++){

                IloEnv env;

                try {

                    instanceReader(0,vert, gclass,ctype, intances_types[i]);
                    
                    // create ILP problem
                    IloModel model(env);
                    cout << "[INFO]: Create variables" << endl;
                    BoolVar3Matrix x(env,num_vertices);
                    BoolVar3Matrix y(env,num_vertices);
                    BoolVar4Matrix f(env,num_vertices);
                    IntVarMatrix lambda(env,num_vertices);
                    allocVars(env, x, y, f, lambda);
                    cout << "[INFO]: Create model" << endl;
                    createModel(model, x, y, f, lambda);
                    IloCplex cplex;
                    cplex = IloCplex(model);
                    //         cplex.exportModel("FEGS.lp");

                    double cpu0, cpu1;
                    cplex.setParam(IloCplex::TiLim, 7200); // tempo limite 2h
                    cplex.setParam(IloCplex::TreLim, 7000); // limite de 7GB de memoria
                    
                    cpu0 = get_wall_time();
                    if (!cplex.solve()){
                        env.error() << "[INFO]: Failed to optimize ILP." << endl;
                        throw(-1);
                    }
                    cpu1 = get_wall_time();
                    double runTime = cpu1 - cpu0;

                    cout << "Solution status = " << cplex.getStatus()   << endl;
                    cout << "Solution value  = " << cplex.getObjValue() << endl;
                    cout << "cplex time: " << cplex.getTime() << endl;
                    cout << "run time: " << runTime << endl;

                    //                displaySolution(cplex,x,y,f,lambda);
                    saveSolution(cplex,x,y,f,ctype);
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
allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda){

    // var x[u][j][s]
    for (int u=0; u<num_vertices; u++){
        x[u] = BoolVarMatrix(env, num_teams);
        for (int j = 0; j < num_teams; j++)
            x[u][j] = IloBoolVarArray(env, num_skills);
    }

    // var y[u][v][j]
    for (int u=0; u<num_vertices; u++){
        y[u] = BoolVarMatrix(env, num_vertices);
        for (int v = 0; v < num_vertices; v++)
            y[u][v] = IloBoolVarArray(env, num_teams);
    }

    // var f[u][v][p][q]
    for (int u=0; u<num_vertices; u++){
        f[u] = BoolVar3Matrix(env, num_vertices);
        for (int v=u+1; v<num_vertices; v++){ //u!=v, u<v
            f[u][v] = BoolVarMatrix(env, num_vertices);
            for (int p=0; p<num_vertices; p++)
                f[u][v][p] = IloBoolVarArray(env, num_vertices);
        }
    }

    // var lambda[u][v]
    for (int u = 0; u < num_vertices; u++){
        lambda[u] = IloIntVarArray(env, num_vertices, 0, num_vertices*(num_vertices-1)); // NECESSARI0 COLOCAR BOUNDS NO INT
    }


}

static void
createModel (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda){

    // add var x
    for(int j=0; j<num_teams; j++) // para todo projeto
        for(int u=0; u< num_vertices; u++) // para todo individuo
            for(int s=0; s<num_skills; s++){ // todas as habilidades que o individuo possui
                if(K[u][s] > 0){ // s(u)
                    char name[30];
                    sprintf(name, "x%d%d%d", u+1, j+1,s+1);
                    x[u][j][s].setName(name);
                    model.add(x[u][j][s]);
                }
            }

    // add var y
    for(int u=0; u<num_vertices; u++) // para todo individuo u
        for(int v = u+1; v<num_vertices; v++) // para todo individuo v, u!=v ou melhor u<v
            for(int j = 0; j<num_teams; j++){ // para todo projeto j
                char name[30];
                sprintf(name, "y%d%d%d", u+1, v+1,j+1);
                y[u][v][j].setName(name);
                model.add(y[u][v][j]);
            }

    // add var f
    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++) // para todo individuo u,v
            if(G[u][v] < 1){ // u,v not in E
                for(int p=0; p<num_vertices; p++)
                    for(int q=0; q<num_vertices; q++) // para todo ARCO p,q
                        if(p!= q && G[p][q] != 0){ // existir aresta -> 2 arcos
                            char name[30];
                            sprintf(name, "F%d%d%d%d", u+1, v+1,p+1, q+1);
                            f[u][v][p][q].setName(name);
                            model.add(f[u][v][p][q]);
                        }
            }

    // add var lambda
    for(int u = 0; u<num_vertices; u++) // para todo individuo u
        for(int v = u+1; v<num_vertices; v++) // para todo individuo v
            if(G[u][v] < 1){ // nao existir aresta positiva entre u e v
                char name[20];
                sprintf(name, "lambd%d%d", u+1, v+1);
                lambda[u][v].setName(name);
                model.add(lambda[u][v]);
            }


    objFunction(model,y,f);
    rest1(model, x);
    rest2(model, x);
    rest3(model,x,y);
    rest4(model,y,f);
    rest5(model,y,f);
    rest6(model,f,lambda);


}

static void
objFunction (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    IloExpr objExpr(env);
    for(int u=0; u<num_vertices; u++)
        for(int v=0; v<num_vertices; v++)
            for(int p=0; p<num_vertices; p++)
                for(int q=0; q<num_vertices; q++)
                    if(G[u][v] < 1 && p!= q && G[p][q] != 0 && u<v){ // cost same team: (u,v) not in E+
                       objExpr += f[u][v][p][q];
                    }
    for(int u=0; u<num_vertices; u++)
        for(int v=0; v<num_vertices; v++)
            for(int j = 0; j<num_teams; j++){
                if(u<v && G[u][v] >= 1){objExpr += y[u][v][j];}  // cost same team: (u,v) in E+
            }

    IloObjective obj = IloMinimize(env, objExpr);
    model.add(obj);
    objExpr.end();
}

static void
rest1 (IloModel model, BoolVar3Matrix x){

    IloEnv env = model.getEnv();

    // each individual use one skill in  maximum
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
            if (R[j][s] > 0) {
                for (int u=0; u < num_vertices; u++) {
                    if(K[u][s] > 0)
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
        for (int v = u+1; v < num_vertices; v++) // u!=v
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
rest4 (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v<num_vertices; v++) // for all u!=v indiduals such that (u,v) not in E+
            if(G[u][v] < 1){
                for(int q = 0; q<num_vertices; q++){ // for fixed vertice q
                    IloExpr sumArcs(env);
                    for(int p = 0; p<num_vertices; p++)// for every arc (p,q) <=> (p,q) in E <=> (q,p) in E
                        if(p!=q && G[p][q] != 0){ //  (p,q) in E <=> (q,p) in E
                            sumArcs += f[u][v][p][q];
                            sumArcs += -f[u][v][q][p];
                        }
                    if(q == u){
                        IloExpr expr(env);
                        for (int j = 0; j < num_teams; j++) {
                            expr += y[u][v][j];
                        }
                        model.add(sumArcs == -expr); expr.end();
                    }else if(q == v){
                        IloExpr expr(env);
                        for (int j = 0; j<num_teams; j++) {
                            expr += y[u][v][j];
                        }
                        model.add(sumArcs == expr); expr.end();
                    }else{
                        model.add(sumArcs == 0);
                    }
                }
            }
}

static void
rest5 (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v<num_vertices; v++)
            if(G[u][v] < 1){
                for(int p=0; p<num_vertices; p++)
                    for(int q=0; q<num_vertices; q++)
                        if(p!=q && G[p][q] != 0){ // (p,q) in E
                            IloExpr expr(env);
                            for (int j=0; j <num_teams; j++) {
                                expr += y[u][v][j];
                            }
                            model.add(f[u][v][p][q] <= expr);
                            expr.end();
                        }
        }
}

static void
rest6 (IloModel model, BoolVar4Matrix f, IntVarMatrix lambda){
    IloEnv env = model.getEnv();

    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v< num_vertices; v++){
            if(G[u][v] < 1 ){ // (u,v) not in E+
                IloExpr expr(env);
                for(int p = 0; p<num_vertices; p++)
                    for(int q = 0; q<num_vertices; q++)
                        if(p!=q && G[p][q] < 0){ // sum negatives arcs(A-)
                            expr += f[u][v][p][q];
                        }
                model.add(expr - 2*lambda[u][v] == 0); // 2k, k in N
                expr.end();
            }
        }
}


static void
displaySolution(IloCplex& cplex,
               BoolVar3Matrix x,
               BoolVar3Matrix y,
               BoolVar4Matrix f,
               IntVarMatrix lambda){



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
    }/*else{
        printf("\n abriu arqv");
    }*/

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
    sprintf(arq, "%s/results/%d_Vertices.ods",CURRENT_DIR, num_vertices);
    
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





