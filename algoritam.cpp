#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <gmp.h>
#include <GLFW/glfw3.h>

using namespace std;
using namespace std::chrono;


struct Graph {
    int U;
    int V;
    int E;
    vector<vector<pair<int,bool>>> edges;

    Graph(int U=0,int V=0,int E=0) : U(U),V(V),E(E) {
        edges.resize(U);
    }

    void addEdge(int u, int v) {
        edges[u].push_back({v,false});
    }
};

vector<vector<int>> generateRandomWeights(Graph G){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1,2*G.E);

    vector<vector<int>> weights(G.U,vector<int>(G.V,0));

    for(int u=0;u<G.U;u++){
        for(pair<int,bool>& edge:G.edges[u]){
            int v=edge.first;
            weights[u][v]=dis(gen);
        }
    }

    return weights;
}

Graph readBipartiteGraphFromFile(const string& filename){
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Unable to open file" << endl;
        exit(EXIT_FAILURE);
    }

    int numVerticesA, numVerticesB, numEdges;
    file >> numVerticesA >> numVerticesB >> numEdges;

    Graph G(numVerticesA, numVerticesB,numEdges);

    for (int i = 0; i < G.E; ++i) {
        int u, v;
        file >> u >> v;
        G.addEdge(u, v);
    }   

    file.close();

    return G;
}

Graph G;
int matchingIterations = 0;
double totalTime = 0.0;

Graph readBipartiteGraphFromStdin() {
    int numVerticesA, numVerticesB, numEdges;
    std::cout << "Unesite broj čvorova u prvoj i drugoj partiji i broj ivica (format: U V E):\n";
    cin >> numVerticesA >> numVerticesB >> numEdges;

    Graph G(numVerticesA, numVerticesB, numEdges);

    std::cout << "Unesite ivice grafa (format: U V):\n";
    for (int i = 0; i < G.E; ++i) {
        int u, v;
        cin >> u >> v;
        G.addEdge(u, v);
    }

    return G;
}

void determinant(mpz_t** mat, int size, mpz_t det){
    if (size == 1) {
        mpz_set(det, mat[0][0]);
    } else if (size == 2) {
        mpz_mul(det, mat[0][0], mat[1][1]);
        mpz_t temp;
        mpz_init(temp);
        mpz_mul(temp, mat[0][1], mat[1][0]);
        mpz_sub(det, det, temp);
        mpz_clear(temp);
    } else {
        mpz_t temp;
        mpz_init(temp);
        mpz_init(det);
        mpz_set_si(det, 0);

        mpz_t sign;
        mpz_init(sign);
        mpz_set_si(sign, 1);

        mpz_t** submatrix = new mpz_t*[size - 1];
        for (int i = 0; i < size - 1; ++i) {
            submatrix[i] = new mpz_t[size - 1];
            for (int j = 0; j < size - 1; ++j) {
                mpz_init(submatrix[i][j]);
            }
        }

        for(int i=0;i<size;i++){
            int subi=0;

            for(int j=1;j<size;j++){
                int subj=0;

                for(int k=0;k<size;k++){
                    if(k!=i){
                        mpz_set(submatrix[subi][subj],mat[j][k]);
                        subj++;
                    }
                }
                subi++;
            }
            determinant((mpz_t**)submatrix, size - 1, temp);
            mpz_mul(temp, mat[0][i], temp);
            mpz_mul(temp, sign, temp);
            mpz_add(det, det, temp);
            mpz_mul_si(sign, sign, -1);
        }
        mpz_clear(temp);
        mpz_clear(sign);

        for (int i = 0; i < size - 1; ++i) {
            for (int j = 0; j < size - 1; ++j) {
                mpz_clear(submatrix[i][j]);
            }
            delete[] submatrix[i];
        }
        delete[] submatrix;
    }
}

void adjointMatrix(mpz_t** mat, int size, mpz_t** adj) {
    if (size <= 1) {
        return;
    }

    mpz_t*** submatrix = new mpz_t**[size];
    for (int i = 0; i < size; ++i) {
        submatrix[i] = new mpz_t*[size - 1];
        for (int j = 0; j < size - 1; ++j) {
            submatrix[i][j] = new mpz_t[size - 1];
            for (int k = 0; k < size - 1; ++k) {
                mpz_init(submatrix[i][j][k]);
            }
        }
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            mpz_t sign, tempDet;
            mpz_init(sign);
            mpz_init(tempDet);

            mpz_set_si(sign, (i + j) % 2 == 0 ? 1 : -1);

            int subi = 0;
            for (int x = 0; x < size; ++x) {
                if (x != i) {
                    int subj = 0;
                    for (int y = 0; y < size; ++y) {
                        if (y != j) {
                            mpz_set(submatrix[i][subi][subj], mat[x][y]);
                            subj++;
                        }
                    }
                    subi++;
                }
            }

            determinant((mpz_t**)submatrix[i], size - 1, tempDet);

            mpz_mul(tempDet, tempDet, sign);
            mpz_set(adj[j][i], tempDet);

            mpz_clear(sign);
            mpz_clear(tempDet);
        }
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size - 1; ++j) {
            for (int k = 0; k < size - 1; ++k) {
                mpz_clear(submatrix[i][j][k]);
            }
            delete[] submatrix[i][j];
        }
        delete[] submatrix[i];
    }
    delete[] submatrix;
}

void perfectMatching(Graph& G, mpz_t** B, int size) {

    mpz_t determinantB;
    mpz_init(determinantB);
    determinant(B, size, determinantB);

    int matchingWeight = 0;
    mpz_t det;
    mpz_init(det);
    mpz_abs(det, determinantB);

    while (mpz_even_p(det) != 0) {
        matchingWeight++;
        mpz_divexact_ui(det, det, 2);
    }
    mpz_t** adjMatrix = new mpz_t*[size];
    for (int i = 0; i < size; ++i) {
        adjMatrix[i] = new mpz_t[size];
        for(int j=0;j<size;j++){
            mpz_init(adjMatrix[i][j]);
        }
    }

    adjointMatrix(B, size, adjMatrix);


    for(int u=0;u<G.U;u++){
       for(pair<int,bool>& edge:G.edges[u]){
            int v=edge.first;
            mpz_t val;
            mpz_init(val);
            mpz_mul(val, adjMatrix[v][u], B[u][v]);
            mpz_tdiv_q_2exp(val, val, matchingWeight);
            if (mpz_odd_p(val) != 0) {
                edge.second = true;
            } else {
                edge.second = false;
            }

            mpz_clear(val);
        }
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            mpz_clear(adjMatrix[i][j]);
        }
        delete[] adjMatrix[i];
    }
    delete[] adjMatrix;


}

bool calculateMatching = false;

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        calculateMatching = true;
    }
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window,GL_TRUE);
    }
}


void display(GLFWwindow* window,mpz_t** B) {
    int width,height;
    glfwGetFramebufferSize(window,&width,&height);

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(1.0,1.0,1.0,1.0);
    glLoadIdentity();


    double rastA = 2.0 / G.U;
    double rastB = 2.0 / G.V;

    glPointSize(10.0);

    for(int u=0;u<G.U;u++){
        glBegin(GL_POINTS);
        glColor3f(0.0,0.0,0.0);
        glVertex2f(-0.4f,-u*rastA+0.9);
        glEnd();
    }
    for(int v=0;v<G.V;v++){
        glBegin(GL_POINTS);
        glColor3f(0.0,0.0,0.0);
        glVertex2f(0.4f,-v*rastB+0.9);
        glEnd();

    }

    for (int u = 0; u < G.U; ++u) {
        for(pair<int,bool>& edge:G.edges[u]){
            int v=edge.first;
            if(edge.second==false){
                glLineWidth(2.0);
                glColor3f(0.0,0.0,0.0);
                glBegin(GL_LINES);
                glVertex2f(-0.4f, -u*rastA+0.9); 
                glVertex2f(0.4f,  -v*rastB+0.9); 
                glEnd();

            }
            else{
                glLineWidth(5.0);
                glColor3f(1.0,0.0,0.0);
                glBegin(GL_LINES);
                glVertex2f(-0.4f, -u*rastA+0.9); 
                glVertex2f(0.4f,  -v*rastB+0.9); 
                glEnd();

            }       
        }
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
}

int main(int argc,char** argv){

    if(argc==2){
        G=readBipartiteGraphFromFile(argv[1]);
    }
    else{
        G=readBipartiteGraphFromStdin();
    }

    GLFWwindow* window;

    if(!glfwInit()){
        return -1;
    }

    window=glfwCreateWindow(800,600,"Bipartitni graf",NULL,NULL);
    if(!window){
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    glfwSetKeyCallback(window, keyCallback);


    vector<vector<int>> D(G.U,vector<int>(G.V,0));
    for(int u=0;u<G.U;u++){
       for(pair<int,bool>& edge:G.edges[u]){
            int v=edge.first;
            D[u][v]=1;
        }
    }

    vector<vector<int>> weights;

    mpz_t** B=new mpz_t*[G.U];
    for (int i = 0; i < G.U; ++i) {
        B[i] = new mpz_t[G.V];
        for (int j = 0; j < G.V; ++j) {
            mpz_init(B[i][j]);
        }
    }

    int matchingCounter=0;
    auto totalTime=0;

    while (!glfwWindowShouldClose(window))
    {

        display(window,B);

        if (calculateMatching) {
            auto start = high_resolution_clock::now();
            matchingCounter++;

            weights = generateRandomWeights(G);

            mpz_t factor;
            mpz_init(factor);
            mpz_set_ui(factor, 2);

            for (int u = 0; u < D.size(); u++) {
                for (int v = 0; v < D[u].size(); v++) {
                    mpz_set_si(B[u][v], D[u][v]);
                    mpz_mul_2exp(B[u][v], B[u][v], weights[u][v]);
                }
            }

            mpz_clear(factor);

            perfectMatching(G, B, G.U);

            calculateMatching = false;
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            totalTime+=duration.count();
            cout << "Matching calculated in " << duration.count() << " microseconds.\n";

        }

        glfwPollEvents();
    }

    cout<<"Total matching count iterations: "<<matchingCounter<<"\n";
    cout<<"Total time: "<<totalTime<<" microseconds\n";
    cout<<"Average time: "<<totalTime/matchingCounter<<" microseconds\n";

    for (int i = 0; i < G.U; ++i) {
        for (int j = 0; j < G.V; ++j) {
            mpz_clear(B[i][j]);
        }
        delete[] B[i];
    }
    delete[] B; 

    glfwTerminate();

    return 0;
}