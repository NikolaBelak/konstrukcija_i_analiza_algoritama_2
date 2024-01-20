#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <GLFW/glfw3.h>

using namespace std;

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

long long int determinant(const vector<vector<int>>& mat){
    int n = mat.size();
    if(n==1){
        return mat[0][0];
    }
    else if(n==2){
        return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
    }

    long long int det=0;
    vector<vector<int>> submatrix(n-1,vector<int>(n-1,0));

    for(int i=0;i<n;i++){
        int subi=0;

        for(int j=1;j<n;j++){
            int subj=0;

            for(int k=0;k<n;k++){
                if(k!=i){
                    submatrix[subi][subj]=mat[j][k];
                    subj++;
                }
            }
            subi++;
        }

        int sign=(i%2==0) ? 1 : -1;
        det+=sign*mat[0][i]*determinant(submatrix);
    }
    return det;
}

vector<vector<long long int>> adjointMatrix(vector<vector<int>>& mat){
    int n=mat.size();
    vector<vector<long long int>> adj(n,vector<long long int>(n,0));

    if(n<=1){
        return adj;
    }
    vector<vector<int>> submatrix(n-1,vector<int>(n-1,0));

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            int subi=0;
            for(int x=0;x<n;x++){
                if(x!=i){
                    int subj=0;
                    for(int y=0;y<n;y++){
                        if(y!=j){
                            submatrix[subi][subj]=mat[x][y];
                            subj++;
                        }
                    }
                    subi++;
                }
            }
            int sign=((i+j)%2==0) ? 1 : -1;
            adj[j][i]=sign*determinant(submatrix);
        }
    }
    return adj;
}

void perfectMatching(Graph& G, vector<vector<int>>& B) {
    long long int determinantB = determinant(B);
    int matchingWeight=0;
    long long int det=abs(determinantB);
    if(det!=0){
        while(det%2==0){
            matchingWeight++;
            det/=2;
        }
    }

    vector<vector<long long int>> adjMatrix=adjointMatrix(B);

    for(int u=0;u<G.U;u++){
       for(pair<int,bool>& edge:G.edges[u]){
            int v=edge.first;
            int val=adjMatrix[v][u]*B[u][v]/pow(2,matchingWeight);
            if(val%2!=0){
                edge.second=true;
            }
        }
    }
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


void display(GLFWwindow* window,vector<vector<int>> B) {
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

    vector<vector<int>> D(G.U,vector<int>(G.V,0));
    for(int u=0;u<G.U;u++){
       for(pair<int,bool>& edge:G.edges[u]){
            int v=edge.first;
            D[u][v]=1;
        }
    }

    vector<vector<int>> weights;
    vector<vector<int>> B(G.U, vector<int>(G.V, 0));

    glfwSetKeyCallback(window, keyCallback);

    while (!glfwWindowShouldClose(window))
    {
        display(window,B);

        if (calculateMatching) {

            weights = generateRandomWeights(G);

            for (int u = 0; u < D.size(); u++) {
                for (int v = 0; v < D[u].size(); v++) {
                    B[u][v] = D[u][v] * pow(2, weights[u][v]);
                }
            }

            perfectMatching(G, B);

            calculateMatching = false;
        }


        glfwPollEvents();
    }

    glfwTerminate();

    return 0;
}