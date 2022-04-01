#include<string>

using namespace std;

class stack
{
    int top,size;
    string arr[100];

public:
    stack():top(-1),size(100){}

    int isfull();
    int isempty();
    void show();
    void push(string data);
    string pop();
    string peek();
};

class basicfns
{
public:
    float trig(int choice,float angdeg);
    float inversetrig(int choice,float val);
    float logfns(int choice,float val);
};

class expeval:public stack,public basicfns
{
    string inexp,postexp[100];
    int outind;
    float result;

    int isoperator(char c);
    int prec(char c);

    float expressioneval(string postexp[],int cnt);

public:
    expeval():outind(-1){};

    void dispinst();
    void intopost(string inexp);
    void run();
};

class quadroot;

class comp
{
    float real;
    float im;

    float mod();
    comp conjugate();

public:
    comp(float r=0,float i=0):real(r),im(i){}

    comp operator+(comp &c);
    comp operator-(comp &c);
    comp operator*(comp &c);
    comp operator/(comp &c);

    friend istream& operator>>(istream &is,comp &c);
    friend ostream& operator<<(ostream &os,comp &c);
    friend int compmenu(comp &c1,comp &c2);
    friend quadroot;
};

class quadroot
{
    int a,b,c,type;
    float x1,x2;
    comp c1,c2;

    public:
    void compute();
    void disppoly();

    friend istream& operator>>(istream &is,quadroot &q);
    friend ostream& operator<<(ostream &os,quadroot &q);
};

class reg
{
protected:
    double x[100];
    double y[100];
    double ypred[500];
    int num;
    string exp;
public:
    reg():exp(""){};
    void input();
    void disp();
    void plot();
    virtual void calcmata()=0;
    virtual void calcmatb()=0;
    virtual void calcparam()=0;
    virtual void calcpredy()=0;
};

class linreg:public reg
{
    float mata[5][5],matb[5];
    float b[5];
public:
    virtual void calcmata();
    virtual void calcmatb();
    virtual void calcparam();
    virtual void calcpredy();
    friend float sum(float arr[],int n,int pow);
};

class quadreg:public reg
{
    float mata[5][5],matb[5];
    float b[5];

public:
    virtual void calcmata();
    virtual void calcmatb();
    virtual void calcparam();
    virtual void calcpredy();
    friend float sum(float arr[],int n,int pow);
};

class cubreg:public reg
{
    float mata[5][5],matb[5];
    float b[5];

public:
    virtual void calcmata();
    virtual void calcmatb();
    virtual void calcparam();
    virtual void calcpredy();
    friend float sum(float arr[],int n,int pow);
};

