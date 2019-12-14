#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <malloc.h>
#include </usr/local/include/eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;
double myPow(double x, int b)
{
    double ans = 1;
    for (int i = 0; i < b; i++)
    {
        ans = ans * x;
    }
    return ans;
}
bool isEqualForVector(vector<int> a, vector<int> b)
{
    if (a.size() != b.size())
    {
        cout << "Error!Two vectors don't have the same size!" << endl;
        exit(1);
    }
    else
    {
        int length = a.size();
        for (int i = 0; i < length; i++)
        {
            if (a[i] != b[i])
            {
                return false;
            }
        }
        return true;
    }
}
class MyMatrix
{
private:
    int xLength, yLength;
    double detEigen()
    {
        MatrixXd tempEigen(this->getxLength(), this->getyLength());
        for (int i = 0; i < this->getxLength(); i++)
        {
            for (int j = 0; j < this->getyLength(); j++)
            {
                tempEigen(i, j) = this->matrix[i][j];
            }
        }
        double deter = tempEigen.determinant();
        return deter;
    }
    double olddet(double **D, int n)
    {
        /*
        This method is possibly outdated since it's too slow.
        */
        cout << "n:" << n << endl;
        /*
        Ref:double **D,int n
        Not modified,just fro ref.
        */
        double d = 0;

        // 一阶二阶直接计算
        if (n == 1)
            d = D[0][0];
        if (n == 2)
            d = D[0][0] * D[1][1] - D[0][1] * D[1][0];
        else
        {
            for (int k = 0; k < n; k++)
            {
                // 为代数余子式申请内存
                double **M;
                M = (double **)malloc((n - 1) * sizeof(double *));
                for (int i = 0; i < n - 1; i++)
                    M[i] = (double *)malloc((n - 1) * sizeof(double));

                // 为代数余子式赋值
                for (int i = 0; i < n - 1; i++)
                    for (int j = 0; j < n - 1; j++)
                        M[i][j] = D[i + 1][j < k ? j : j + 1];

                // 按第一行展开，递归计算行列式，注意元素0则不展开可以加快计算速度
                if (D[0][k])
                    d += D[0][k] * olddet(M, n - 1) * (((2 + k) % 2) ? -1 : 1);

                // 释放内存
                for (int i = 0; i < n - 1; i++)
                    free(M[i]);
                free(M);
            }
        }
        return d;
    }
    MyMatrix getAccompany()
    {
        MatrixXd tempEigen(this->getxLength(), this->getyLength());
        for (int i = 0; i < this->getxLength(); i++)
        {
            for (int j = 0; j < this->getyLength(); j++)
            {
                tempEigen(i, j) = this->matrix[i][j];
            }
        }
        MatrixXd accompanyAns = tempEigen.adjoint();
        MyMatrix temp(xLength, yLength, "zeros");
        for (int i = 0; i < this->getxLength(); i++)
        {
            for (int j = 0; j < this->getyLength(); j++)
            {
                temp.matrix[i][j] = accompanyAns(i, j);
            }
        }
        return temp;
    }
    MyMatrix oldgetAccompany()
    {
        /*
        Outdated
        */
        MyMatrix ans(xLength, yLength, "zeros");
        if (xLength == 1)
        {
            ans.matrix[0][0] = 1;
            return ans;
        }
        int i, j, k, t;
        int n = xLength;
        MyMatrix temp(xLength - 1, yLength - 1, "zeros");
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                for (k = 0; k < n - 1; k++)
                {
                    for (t = 0; t < n - 1; t++)
                    {
                        temp.matrix[k][t] = this->matrix[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
                    }
                }
                ans.matrix[j][i] = temp.getDet();
                if ((i + j) % 2 == 1)
                {
                    ans.matrix[j][i] = -ans.matrix[j][i];
                }
            }
        }
        return ans;
    }

public:
    double **matrix = nullptr;
    MyMatrix(int xlength, int ylength, string type)
    {
        xLength = xlength;
        yLength = ylength;
        matrix = new double *[xLength];
        for (int i = 0; i < xLength; i++)
        {
            matrix[i] = new double[yLength];
        }
        if (type == "zeros")
        {
            for (int i = 0; i < xLength; i++)
            {
                for (int j = 0; j < yLength; j++)
                {
                    matrix[i][j] = 0;
                }
            }
        }
        else if (type == "ones")
        {
            for (int i = 0; i < xLength; i++)
            {
                for (int j = 0; j < yLength; j++)
                {
                    if (i == j)
                        matrix[i][j] = 1;
                    else
                    {
                        matrix[i][j] = 0;
                    }
                }
            }
        }
        else
        {
            cout << "Error:the type of matrix to be initialized is undefined!" << endl;
        }
    }
    int getxLength()
    {
        return xLength;
    }
    int getyLength()
    {
        return yLength;
    }
    ~MyMatrix()
    {
        for (int i = 0; i < xLength; i++)
        {
            delete[] matrix[i];
        }
        delete[] matrix;
    }
    MyMatrix operator+(const MyMatrix &b)
    {
        MyMatrix temp(this->getxLength(), this->getyLength(), "zeros");
        if (b.xLength == this->getxLength() && b.yLength == this->getyLength())
        {
            for (int i = 0; i < this->getxLength(); i++)
            {
                for (int j = 0; j < this->getyLength(); j++)
                {
                    temp.matrix[i][j] = this->matrix[i][j] + b.matrix[i][j];
                }
            }
            return temp;
        }
        else
        {
            cout << "the lengths of two matrixes in dimensions are different!" << endl;
            exit(1);
        }
    }
    MyMatrix operator-(const MyMatrix &b)
    {
        MyMatrix temp(this->getxLength(), this->getyLength(), "zeros");
        if (b.xLength == this->getxLength() && b.yLength == this->getyLength())
        {
            for (int i = 0; i < this->getxLength(); i++)
            {
                for (int j = 0; j < this->getyLength(); j++)
                {
                    temp.matrix[i][j] = this->matrix[i][j] - b.matrix[i][j];
                }
            }
            return temp;
        }
        else
        {
            cout << "Error:the lengths of two matrixes in dimensions are different!" << endl;
            exit(1);
        }
    }
    MyMatrix operator*(const MyMatrix &b)
    {
        MyMatrix temp(this->getxLength(), b.yLength, "zeros");
        if (b.xLength == this->getyLength())
        {
            for (int i = 0; i < this->getxLength(); i++)
            {
                for (int j = 0; j < b.yLength; j++)
                {
                    for (int k = 0; k < this->getyLength(); k++)
                        temp.matrix[i][j] += (double)this->matrix[i][k] * (double)b.matrix[k][j];
                }
            }
            return temp;
        }
        else
        {
            cout << "Error:the first matrix's yLength must be same as the second matrix's xLength!" << endl;
            exit(1);
        }
    }
    MyMatrix operator*(const double &b)
    {
        MyMatrix temp(this->getxLength(), this->getyLength(), "zeros");
        for (int i = 0; i < this->getxLength(); i++)
        {
            for (int j = 0; j < this->getyLength(); j++)
            {
                    temp.matrix[i][j] = this->matrix[i][j] * b;
            }
        }

        return temp;
    }
    MyMatrix operator/(const double b)
    {
        MyMatrix ans(xLength, yLength, "zeros");
        for (int i = 0; i < xLength; i++)
        {
            for (int j = 0; j < yLength; j++)
            {
                ans.matrix[i][j] = this->matrix[i][j] / b;
            }
        }
        return ans;
    }
    void operator=(const MyMatrix &b)
    {
        setEquals(b);
    }
    void setEquals(const MyMatrix &b)
    {
        if (this->getxLength() == b.xLength && this->getyLength() == b.yLength)
        {
            for (int i = 0; i < xLength; i++)
            {
                for (int j = 0; j < yLength; j++)
                {
                    this->matrix[i][j] = b.matrix[i][j];
                }
            }
            return;
        }
        else
        {
            cout << "error!:Two matrixes must have the same length in every dimension" << endl;
            exit(1);
        }
    }
    MyMatrix getConvertMyMatrix()
    {
        MyMatrix temp(this->getyLength(), this->getxLength(), "zeros");
        for (int i = 0; i < this->getxLength(); i++)
        {
            for (int j = 0; j < this->getyLength(); j++)
            {
                temp.matrix[j][i] = this->matrix[i][j];
            }
        }
        return temp;
    }
    void showMyMatrix()
    {
        for (int i = 0; i < this->getxLength(); i++)
        {
            for (int j = 0; j < this->getyLength(); j++)
            {
                cout << matrix[i][j] << '\t';
            }
            cout << endl;
        }
    }
    double getDet()
    {
        /*
        Ref:https://blog.csdn.net/Eyizoha/article/details/89376301
        Modified.
        */
        if (xLength != yLength)
        {
            cout << "Error!xLength must be equal to yLength to calculate det" << endl;
            exit(1);
        }
        double ans = detEigen();
        return ans;
    }
    MyMatrix getInverse()
    {
        MyMatrix inverse = this->getAccompany();
        inverse = inverse / this->getDet();
        return inverse;
    }
    MyMatrix subMyMatrixByRow(int row)
    {
        MyMatrix ans(1, this->getyLength(), "zeros");
        for (int j = 0; j < this->getyLength(); j++)
        {
            ans.matrix[0][j] = this->matrix[row][j];
        }
        return ans;
    }
};
class ReadData
{
private:
    vector<double> XList;
    vector<double> YList;

public:
    void readFile(string file)
    {
        ifstream in(file);
        if (!in.is_open())
        {
            cout << "error!The file can't be opened!" << endl;
            exit(1);
        }
        double temp;
        while (!in.eof())
        {
            in >> temp;
            in >> temp;
            XList.push_back(temp);
            in >> temp;
            YList.push_back(temp);
        }
        in.close();
        XList.pop_back();
        YList.pop_back();
        return;
    }
    int getLength()
    {
        return XList.size();
    }
    vector<double> getXList()
    {
        return XList;
    }
    vector<double> getYList()
    {
        return YList;
    }
    MyMatrix getDataInMyMatrix()
    {
        MyMatrix ans(this->getLength(), 2, "zeros");
        for (int i = 0; i < this->getLength(); i++)
        {
            ans.matrix[i][0] = XList[i];
            ans.matrix[i][1] = YList[i];
        }
        return ans;
    }
};
class GMM
{
private:
    int modelNumber;
    int countExecuteTime = 0;
    vector<MyMatrix *> ConverianceList;
    double *piror;
    MyMatrix *average;
    double *weight;
    MyMatrix *componentPossiblity;
    vector<int> classify()
    {
        vector<int> ansSet;
        for (int i = 0; i < componentPossiblity->getxLength(); i++)
        {
            double maxmium = 0;
            int index = 0;
            for (int j = 0; j < modelNumber; j++)
            {
                if (componentPossiblity->matrix[i][j] > maxmium)
                {
                    maxmium = componentPossiblity->matrix[i][j];
                    index = j;
                }
            }
            ansSet.push_back(index);
        }
        return ansSet;
    }
    double calculateN(MyMatrix &data, int i, int k)
    {
        double ans = 0;
        MyMatrix temp(1, data.getyLength(), "zeros");
        for (int z = 0; z < data.getyLength(); z++)
        {
            temp.matrix[0][z] = data.matrix[i][z];
        }
        /*
        cout << "data:" << endl;
        temp.showMyMatrix();
        */
        temp = temp - average->subMyMatrixByRow(k);
        /*
        cout << "temp:" << endl;
        temp.showMyMatrix();
        */
        MyMatrix temp2 = temp.getConvertMyMatrix();
        MyMatrix *tempConveriance = ConverianceList[k];
        /*
        cout << "tempConveriance:" << endl;
        tempConveriance->showMyMatrix();
        */
        /*
        have to take the equation apart to find the bug
        */
        double ansP1 = myPow(2 * M_PI, modelNumber);
        ansP1 = ansP1 * tempConveriance->getDet();
        ansP1 = 1 / sqrt(ansP1);
        /*
        cout << "tempConveriance:" << endl;
        tempConveriance->showMyMatrix();
        */
        MyMatrix ansP2 = tempConveriance->getInverse() * temp2;
        /*
        cout << "ansP2:" << endl;
        ansP2.showMyMatrix();
        cout << "temp:" << endl;
        temp.showMyMatrix();
        */
        MyMatrix ansP3 = temp * ansP2;
        /*
        cout << "ansP3:" << endl;
        ansP3.showMyMatrix();
        */
        double ansP4 = ansP3.getDet();
        ans = (ansP1 * exp(-0.5 * ansP4));
       // cout<<"the answer of calculating N is:"<<ans<<endl;
        return ans;
    }
    void Epart(MyMatrix &data)
    {
        countExecuteTime++;
        /*
        cout << "execute" << countExecuteTime << "times\nweight=" << weight[0] << "," << weight[1] << endl;
        cout << "average:" << endl;
        average->showMyMatrix();
        cout << "converiance:" << endl;
        ConverianceList[0]->showMyMatrix();
        ConverianceList[1]->showMyMatrix();*/
        /*cout << "componentPossiblity(before change):" << endl;
        componentPossiblity->showMyMatrix();*/
        /*
        problem exists in the epart
        */
        for (int i = 0; i < data.getxLength(); i++)
        {
            /*if (i == data.getxLength() - 1)
            {
                cout << "the last row:" << endl;
            }*/
            double toDivision = 0;
            for (int k = 0; k < modelNumber; k++)
            {
                componentPossiblity->matrix[i][k] = weight[k] * calculateN(data, i, k);
                /*cout << k + 1 << "-th is:" << endl;
                cout << calculateN(data, i, k) << endl;*/
                toDivision += weight[k] * calculateN(data, i, k);
            }
            for (int k = 0; k < modelNumber; k++)
            {
                if (toDivision == 0)
                {
                    if (componentPossiblity->matrix[i][k] != 0)
                    {
                        cout<<"error:zero can't be divided"<<endl;
                        exit(1);
                    }
                }
                else{
                    componentPossiblity->matrix[i][k] = componentPossiblity->matrix[i][k] / toDivision;
                    if(componentPossiblity->matrix[i][k]<0.0001){
                        componentPossiblity->matrix[i][k]=0;
                    }
                }
                    
            }
        }
        /*
        cout << "componentPossiblity(after change):" << endl;
        componentPossiblity->showMyMatrix();*/
    }
    void Mpart(MyMatrix &data)
    {
        //update the average
        for (int k = 0; k < modelNumber; k++)
        {
            double Nk = 0;
            for (int i = 0; i < data.getxLength(); i++)
            {
                Nk += componentPossiblity->matrix[i][k];
            }
            for (int j = 0; j < data.getyLength(); j++)
            {
                double tempSum = 0;
                for (int i = 0; i < data.getxLength(); i++)
                {
                    tempSum += componentPossiblity->matrix[i][k] * data.matrix[i][j];
                }
                average->matrix[k][j] = (1 / Nk) * tempSum;
            }
            //update the converiance
            MyMatrix forConveriance(data.getyLength(), data.getyLength(), "zeros");
            for (int i = 0; i < data.getxLength(); i++)
            {
                MyMatrix Xi(1, data.getyLength(), "zeros");
                for (int j = 0; j < data.getyLength(); j++)
                {
                    Xi.matrix[0][j] = data.matrix[i][j] - average->matrix[k][j];
                }
                MyMatrix XiTran = Xi.getConvertMyMatrix();
                /*
                cout<<"Xi:"<<endl;
                Xi.showMyMatrix();
                cout<<"Xitran:"<<endl;
                XiTran.showMyMatrix();
                cout<<"poss:"<<componentPossiblity->matrix[i][k]<<endl;
                cout<<"subSum:"<<endl;
                (XiTran * Xi).showMyMatrix();
                cout<<"sum:"<<endl;
                */
                //((XiTran * Xi) * componentPossiblity->matrix[i][k]).showMyMatrix();
                forConveriance = forConveriance + (XiTran * Xi) * componentPossiblity->matrix[i][k];
            }
            
            
            //forConveriance.showMyMatrix();

             MyMatrix *ans = new MyMatrix(data.getyLength(), data.getyLength(), "zeros");
            ans->setEquals((forConveriance / Nk));
            ConverianceList[k] = ans;
            //update the weight
            weight[k] = Nk / data.getxLength();
        }
        return;
    }

public:
    GMM(int MODELNUMBER)
    {
        this->modelNumber = MODELNUMBER;
        piror = new double[MODELNUMBER];
        weight = new double[modelNumber];
    }
    void reportPara(){
        cout<<"-------------------report---------------------------"<<endl;
        cout<<"modelNumber:"<<modelNumber<<endl;
        cout<<"ConverianceList"<<endl;
        for(int i=0;i<modelNumber;i++){
            ConverianceList[i]->showMyMatrix();
        }
        cout<<"piror:"<<endl;
        for(int i=0;i<modelNumber;i++){
            cout<<piror[i];
        }
        cout<<"average:"<<endl;
        average->showMyMatrix();
        cout<<"weight:"<<endl;
                for(int i=0;i<modelNumber;i++){
            cout<<weight[i];
        }
    }
    vector<int> run(MyMatrix &data)
    {
        /*
        the main entry to GMM
        */
        componentPossiblity = new MyMatrix(data.getxLength(), modelNumber, "zeros");
        this->ConverianceList.clear();
        for (int i = 0; i < modelNumber; i++)
        {
            this->ConverianceList.push_back(new MyMatrix(data.getyLength(), data.getyLength(), "ones"));
        }
        this->average = new MyMatrix(modelNumber, data.getyLength(), "zeros");
        srand(time(NULL));
        /*
        update the weight for 2 models*/
        weight[0] = 0.5;
        weight[1] = 0.5;
        for (int i = 0; i < modelNumber; i++)
        {
            //weight[i] = 1.0 / modelNumber;
            int randomNumber = rand() % data.getxLength();
            for (int k = 0; k < data.getyLength(); k++)
            {
                average->matrix[i][k] = data.matrix[randomNumber][k];
            }
        }
        int CYCLETIME = 2;
        for (int i = 0; i < CYCLETIME; i++)
        {
            Epart(data);
            Mpart(data);
        }
        Epart(data);
        Mpart(data);
        vector<int> ans1 = classify();
        Epart(data);
        Mpart(data);
        vector<int> ans2 = classify();
        while (!isEqualForVector(ans1, ans2))
        {
            for (int i = 0; i < CYCLETIME; i++)
            {
                Epart(data);
                Mpart(data);
            }
            Epart(data);
            Mpart(data);
            vector<int> ans1 = classify();
            Epart(data);
            Mpart(data);
            vector<int> ans2 = classify();
        }
        return ans1;
    }
};

int main()
{
    const int MODELNUMBER = 2;
    ReadData rd;
    rd.readFile("/home/jiading/Desktop/C++/GMM/dataGMM.txt");
    MyMatrix data = rd.getDataInMyMatrix();
    GMM mm(MODELNUMBER);
    vector<int> ans = mm.run(data);
    for (int i = 0; i < ans.size(); i++)
    {
        cout << ans[i] << endl;
    }
    mm.reportPara();
}
