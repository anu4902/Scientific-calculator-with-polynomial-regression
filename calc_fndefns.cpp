#include "calc_classdefn.h"

#include "pbPlots.hpp"
#include "supportLib.hpp"
#include<iostream>
#include<iomanip>
#include<string>
#include<stdlib.h>
#include<math.h>
#define PI 3.14
#include <codecvt>
#include <locale>

using namespace std;
using convert_t = std::codecvt_utf8<wchar_t>;
std::wstring_convert<convert_t, wchar_t> strconverter;

std::string to_string(std::wstring wstr)
{
    return strconverter.to_bytes(wstr);
}

std::wstring to_wstring(std::string str)
{
    return strconverter.from_bytes(str);
}

//stack

int stack::isfull()
{
    if(top==size-1)
        return 1;
    else
        return 0;
}

int stack::isempty()
{
    if(top==-1)
        return 1;
    else
        return 0;
}

void stack::show()
{
    int i;
    if(isempty())
        printf("\nStack is empty!\n");

    else
    {
        printf("\nStack: ");
        for(i=top;i>=0;i--)
            cout<<arr[i];
    }
}

void stack::push(string data)
{
    if(isfull())
        printf("Stack is already full!\n");

    else
    {
        arr[++top]=data;
    }
}

string stack::pop()
{
    string a;
    if(isempty())
        printf("\nStack is empty!\n");

    else
    {
        a=arr[top--];
        return a;
    }
}

string stack::peek()
{
    return arr[top];
}

//basicfns

float basicfns::trig(int choice,float angdeg)
{
    float angrad=angdeg*PI/180;

    switch(choice)
    {
    case 1:
        return sin(angrad);

    case 2:
        return cos(angrad);

    case 3:
        return tan(angrad);
    }
}

float basicfns::inversetrig(int choice,float val)
{
    switch(choice)
    {
    case 1:
        return asin(val);
    case 2:
        return acos(val);
    case 3:
        return atan(val);
    }
}

float basicfns::logfns(int choice,float val)
{
    switch(choice)
    {
    case 1:
        return log10(val);
    case 2:
        return log(val);
    }
}



string getString(char x)
{
    string s(1, x);
    return s;
}

//expeval

void expeval::dispinst()
{
    cout<<left;
    cout<<"\tKindly follow the given syntax:\n";
    cout<<"\n->Arithmetic operations:\n\t\t"<<setw(15)<<"Operators"<<": +, -, *, /, ^\n\t\t"<<setw(15)<<"Eg"<<": -3.5*2+2/6^2\n\n";cout<<"->Trigonometric functions (Enter angle measure in degrees):\n\t\t"<<setw(15)<<"Functions"<<": sin, cos, tan, sininv, cosinv, taninv\n\t\t"<<setw(15)<<"Eg"<<": sin(0.5)-2+sininv(0)\n\n";
    cout<<"->Logarithmic functions:\n\t\t"<<setw(15)<<"Functions"<<": log(x)-Base 10 ; ln(x) - Base e\n\t\t"<<setw(15)<<"Eg"<<": -2+log(2.5)-ln(6)\n\n";
    cout<<"->Available constants:\n\t\t"<<setw(15)<<"e = 2.718 "<<"; pi = 3.14\n\t\tEg: 2*pi*3+e^2\n\n";
    cout<<"\t\t(Enter exit() to close the application)\n\n";
}

int expeval::isoperator(char c)
{
   char oplist[6]={'+','-','*','/','^'};
   for(int i=0;i<5;i++)
   {
       if(oplist[i]==c)
            return 1;
   }
   return 0;
}

int expeval::prec(char c)
{
    int pre;
    switch(c)
    {
        case '+':
        case '-':
            pre= 1;
            break;
        case '*':
        case '/':
            pre=2;
            break;
        case '^':
            pre=3;
            break;
    }
    return pre;
}

void expeval::run()
{
    while(1)
    {
        cout<<"\n\n>> ";
        outind=-1;
        cin>>inexp;
        if(inexp!="exit()")
            intopost(inexp);
        else
            break;
    }
}

float expeval::expressioneval(string postexp[],int cnt)
{
    string s;
    float op1,op2,temp;

    for(int i=0;i<=cnt;i++)
    {
        s=postexp[i];

        if(isoperator(s[0]) && s.length()==1)
        {
            op2=stof(pop());
            op1=stof(pop());

            switch(s[0])
            {
                case '+':
                    temp=op1+op2;
                    break;
                case '-':
                    temp=op1-op2;
                    break;
                case '*':
                    temp=op1*op2;
                    break;
                case '/':
                    temp=op1/op2;
                    break;
                case '^':
                    temp=pow(op1,op2);
                    break;
            }
            push(to_string(temp));
        }

        else
        {
            push(s);
        }
    }
    return stof(pop());
}

void expeval::intopost(string inexp)
{
    for(int i=0;i<inexp.length();i++)
    {
        if(isdigit(inexp[i]))
        {
            string s=getString(inexp[i]);
            while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
            {
                s.append(getString(inexp[++i]));
            }
            postexp[++outind]=s;
        }

        else if(isalpha(inexp[i]))
        {
            //sin and sininv
            if(inexp[i]=='s' && inexp[i+1]=='i' && inexp[i+2]=='n')
            {
                if(inexp[i+3]=='(')
                {
                    i+=4;
                    string s=getString(inexp[i]);
                    while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                    {
                        s.append(getString(inexp[++i]));
                    }
                    i++;
                    postexp[++outind]=to_string(trig(1,stof(s)));
                }
                else if(inexp[i+3]=='i'&&inexp[i+4]=='n'&&inexp[i+5]=='v'&&inexp[i+6]=='(')
                {
                    i+=7;
                    string s=getString(inexp[i]);
                    while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                    {
                        s.append(getString(inexp[++i]));
                    }
                    i++;
                    postexp[++outind]=to_string(inversetrig(1,stof(s)));
                }
            }

            //cos and cosinv
            else if(inexp[i]=='c' && inexp[i+1]=='o' && inexp[i+2]=='s')
            {
                if(inexp[i+3]=='(')
                {
                    i+=4;
                    string s=getString(inexp[i]);
                    while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                    {
                        s.append(getString(inexp[++i]));
                    }
                    i++;
                    postexp[++outind]=to_string(trig(2,stof(s)));
                }
                else if(inexp[i+3]=='i'&&inexp[i+4]=='n'&&inexp[i+5]=='v'&&inexp[i+6]=='(')
                {
                    i+=7;
                    string s=getString(inexp[i]);
                    while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                    {
                        s.append(getString(inexp[++i]));
                    }
                    i++;
                    postexp[++outind]=to_string(inversetrig(2,stof(s)));
                }
            }
                //tan and taninv
            else if(inexp[i]=='t' && inexp[i+1]=='a' && inexp[i+2]=='n')
            {
                if(inexp[i+3]=='(')
                {
                    i+=4;
                    string s=getString(inexp[i]);
                    while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                    {
                        s.append(getString(inexp[++i]));
                    }
                    i++;
                    postexp[++outind]=to_string(trig(3,stof(s)));
                }
                else if(inexp[i+3]=='i'&&inexp[i+4]=='n'&&inexp[i+5]=='v'&&inexp[i+6]=='(')
                {
                    i+=7;
                    string s=getString(inexp[i]);
                    while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                    {
                        s.append(getString(inexp[++i]));
                    }
                    i++;
                    postexp[++outind]=to_string(inversetrig(3,stof(s)));
                }
            }

            else if(inexp[i]=='l' && inexp[i+1]=='n' && inexp[i+2]=='(')
            {
                i+=3;
                string s=getString(inexp[i]);
                while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                {
                    s.append(getString(inexp[++i]));
                }
                i++;
                postexp[++outind]=to_string(logfns(2,stof(s)));
            }

            else if(inexp[i]=='l' && inexp[i+1]=='o' && inexp[i+2]=='g' && inexp[i+3]=='(')
            {
                i+=4;
                string s=getString(inexp[i]);
                while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                {
                    s.append(getString(inexp[++i]));
                }
                i++;
                postexp[++outind]=to_string(logfns(1,stof(s)));
            }

            else if(inexp[i]=='e')
            {
               postexp[++outind]="2.718";
            }

            else if(inexp[i]=='p' && inexp[i+1]=='i')
            {
               postexp[++outind]="3.14";
               i++;
            }

            else
            {
                cout<<"Invalid!";
                return;
            }
        }

        else if(inexp[i]=='(')
            push(getString(inexp[i]));

        else if(inexp[i]==')')
        {
            while(peek()!="(")
                  postexp[++outind]=pop();

            pop();
        }

        else if(isoperator(inexp[i]))
        {
            if(inexp[i]=='-' && (i==0||isoperator(inexp[i-1])))
            {
                string s=getString(inexp[i]);
                while(i+1<inexp.length() && (isdigit(inexp[i+1])||inexp[i+1]=='.'))
                {
                    s.append(getString(inexp[++i]));
                }
                postexp[++outind]=s;
            }
            else
            {
                comp:
                if(isempty()||peek()=="(")
                    push(getString(inexp[i]));

                else if(prec(inexp[i])>prec(peek()[0]))
                    push(getString(inexp[i]));
                else
                {
                    postexp[++outind]=pop();
                    goto comp;
                }
            }
        }

        else if(inexp[i]==' ')
            continue;

        else
        {
            cout<<"Invalid!";
            return;
        }
    }

    while(!isempty())
    {
        postexp[++outind]=pop();
    }

    /*cout<<"\n";
    for(int i=0;i<=outind;i++)
        cout<<postexp[i]<<" ";*/

    result=expressioneval(postexp,outind);
    cout<<setprecision(4)<<result;
}

//complex

float comp::mod()
{
    return pow(pow(real,2)+pow(im,2),0.5);
}

comp comp::conjugate()
{
    comp temp;
    temp.real=real;
    temp.im=(-1)*im;
    return temp;
}

comp comp::operator+(comp &c)
{
    comp temp;
    temp.real=real+c.real;
    temp.im=im+c.im;
    return temp;
}

comp comp::operator-(comp &c)
{
    comp temp;
    temp.real=real-c.real;
    temp.im=im-c.im;
    return temp;
}

comp comp::operator*(comp &c)
{
    comp temp;
    temp.real=(real*c.real)-(im*c.im);
    temp.im=(real*c.im)+(im*c.real);
    return temp;
}

comp comp::operator/(comp &c)
{
    comp ct,ct2,temp;
    ct.real=real;
    ct.im=im;
    temp=c.conjugate();
    ct2=ct*temp;
    ct.real=(ct2.real)/(pow(c.mod(),2));
    ct.im=(ct2.im)/(pow(c.mod(),2));
    return ct;
}

istream& operator>>(istream &is,comp &c)
{
    cout<<"Enter real part     : ";
    is>>c.real;
    cout<<"Enter imaginary part: ";
    is>>c.im;
    return is;
}

ostream& operator<<(ostream &os,comp &c)
{
    if(c.real==0 && c.im==0)
        cout<<"0\n";
    else
    {
        if(c.real!=0)
            os<<setprecision(2)<<c.real<<" ";
        if(c.im>0 && c.real!=0)
            cout<<"+ ";
        if(c.im!=0)
            if(c.im==1)
                os<<"i"<<"\n";
            else if(c.im==-1)
                os<<"- i"<<"\n";
            else
                os<<setprecision(2)<<c.im<<"i";
    }
    return os;
}

int compmenu(comp &c1,comp &c2)
{
    int choice,flag=0;
    comp result;
    string res;

    cout<<"\n\n\t\t1. Add\n\t\t2. Subtract\n\t\t3. Multiply\n\t\t4. Divide\n\nEnter choice: ";
    cin>>choice;

    switch(choice)
    {
    case 1:
        result=c1+c2;
        res="Sum";
        break;
    case 2:
        result=c1-c2;
        res="Difference";
        break;
    case 3:
        result=c1*c2;
        res="Product";
        break;
    case 4:
        if(c2.real==0 && c2.im==0)
        {
            cout<<"Cannot divide by zero!";
            flag=1;
            break;
        }
        result=c1/c2;
        res="Quotient";
        break;
    }
    if(flag==0)
    {
        cout<<"\n\nFirst number : "<<c1;
        cout<<"\nSecond number: "<<c2<<"\n";
        cout<<left<<setw(13)<<res<<": "<<result;
    }
    if(choice!=0)
        sleep(1);
    return choice;
}

//quad root finder

void quadroot::disppoly()
{
    cout<<"\nThe quadratic polynomial is: ";
    if(a!=1 && a!=-1)
        cout<<a<<"x^2";
    else if(a==1)
        cout<<"x^2";
    else if(a==-1)
        cout<<"-x^2";

    if(b>0)
        cout<<" + ";
    if(b!=0)
    {
        if(b!=1 && b!=-1)
            cout<<b<<"x";
        else if(b==1)
            cout<<"x";
        else if(b==-1)
            cout<<"-x";
    }

    if(c>0)
        cout<<" + ";
    if(c!=0)
        cout<<c<<"\n";
}

void quadroot::compute()
{
    if(pow(b,2)-(4*a*c)<0)
    {
        type=2;
        c1.real=c2.real=(float)((-1)*b)/(2*a);
        c1.im=pow(abs(pow(b,2)-(4*a*c)),0.5)/(2*a);
        //cout<<"det="<<pow(b,2)-(4*a*c)<<" abs="<<abs(pow(b,2)-(4*a*c))<<" root="<<pow(abs(pow(b,2)-(4*a*c)),0.5);
        c2.im=(-1)*pow(abs(pow(b,2)-(4*a*c)),0.5)/(2*a);
    }

    else
    {
        type=1;
        x1=((-1)*b+pow(pow(b,2)-(4*a*c),0.5))/(2*a);
        x2=((-1)*b-pow(pow(b,2)-(4*a*c),0.5))/(2*a);
    }
}

istream& operator>>(istream &is,quadroot &q)
{
    cout<<"Enter coefficient of x^2 : ";
    is>>q.a;
    while(q.a==0)
    {
        cout<<"\nCoefficient of x^2 cannot be zero!";
        cout<<"Enter coefficient of x^2 : ";
        is>>q.a;
    }
    cout<<"Enter coefficient of x^1 : ";
    is>>q.b;
    cout<<"Enter constant term      : ";
    is>>q.c;
    return is;
}

ostream& operator<<(ostream &os,quadroot &q)
{
    if(q.type==2)
    {
        os<<"\nReal roots doesn't exist.\n\nComplex roots are:\n\t";
        cout<<q.c1<<" and "<<q.c2<<"\n\n";
    }

    else
        os<<"\nRoots are:\n\t\tx1= "<<q.x1<<"  ;   x2= "<<q.x2<<"\n";
    return os;
}

//helper fns
double sum(double arr[],int n,int power)
{
    double sum=0;
    for(int i=0;i<n;i++)
    {
        sum+=pow(arr[i],power);
    }
    return sum;
}

double minarr(double ar[],int n)
{
    double m=ar[0];
    for(int i=0;i<n;i++)
        if(ar[i]<m)
            m=ar[i];
    return m;
}

double maxarr(double ar[],int n)
{
    double m=ar[0];
    for(int i=0;i<n;i++)
        if(ar[i]>m)
            m=ar[i];
    return m;
}

void getCofactor(float mat[5][5], float temp[5][5], int p,int q, int n)
{
    int i = 0, j = 0;

    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            if (row != p && col != q)
            {
                temp[i][j++] = mat[row][col];
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

float determinantOfMatrix(float mat[5][5], int n)
{
    float D = 0;

    if (n == 1)
        return mat[0][0];

    float temp[5][5]; // To store cofactors

    int sign = 1; // To store sign multiplier

    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of mat[0][f]
        getCofactor(mat, temp, 0, f, n);
        D += sign * mat[0][f]
             * determinantOfMatrix(temp, n - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }
    return D;
}

void adjoint(float A[5][5],float adj[5][5])
{
    int N=4;
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1;
    float temp[5][5];

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinantOfMatrix(temp, N-1));
        }
    }
}

bool inverse(float A[5][5], float inverse[5][5])
{
    int N=4;
    // Find determinant of A[][]
    float det = determinantOfMatrix(A, N);
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    float adj[5][5];
    adjoint(A, adj);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            inverse[i][j] = adj[i][j]/float(det);

    return true;
}

//regression
void reg::input()
{
    string s,sx="",sy="";
    cout<<"\nEnter number of points(minimum 5) : ";
    cin>>num;

    while(num<5)
    {
        cout<<"\nEnter number of points(minimum 5) : ";
        cin>>num;
    }
    cout<<"\nEnter points in the format x,y:\n\n";
    for(int i=0,k=0;i<num;i++)
    {
        k=0;
        sx="";
        sy="";
        cout<<"Enter point "<<i+1<<": ";
        cin>>s;
        while(s[k]!=',')
        {
            sx.append(getString(s[k]));
            k++;
        }
        x[i]=stod(sx);
        k++;
        while(k<s.length())
        {
            sy.append(getString(s[k]));
            k++;
        }
        y[i]=stod(sy);
    }
}

void reg::disp()
{
    cout<<left;
    cout<<"\nThe given points are:\n";
    for(int i=0;i<num;i++)
    {
        cout<<"\t("<<x[i]<<" , "<<y[i]<<")\n";
    }
}

void reg::plot()
{
    double xpl[500];
    int cnt;
    //wchar_t=L(
    wstring polystr=to_wstring("y = "+exp);

    cout<<"\nxplot:\n";
    for(double i=minarr(x,num)-5;i<=maxarr(x,num)+5;i+=0.1)
    {
        xpl[cnt]=i;
        //cout<<xpl[cnt]<<endl;
        cnt++;
    }

    int nval=(int)(maxarr(x,num)+5-minarr(x,num)+5)*10;
    vector<double> xvect(xpl,xpl+nval);
    vector<double> yvect(ypred,ypred+nval);

    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
	series->xs = &xvect;
	series->ys = &yvect;
	series->linearInterpolation = true;
	series->lineType = toVector(L"dashed");
	series->lineThickness = 2;
	series->color = GetGray(0.3);

	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 1980;
	settings->height = 1024;
	settings->autoBoundaries = true;
	settings->autoPadding = true;
	settings->title = toVector(polystr.c_str());
	settings->scatterPlotSeries->push_back(series);



    RGBABitmapImageReference *imageref=CreateRGBABitmapImageReference();
    DrawScatterPlotFromSettings(imageref, settings);
    vector<double> *pngdata=ConvertToPNG(imageref->image);
    WriteToFile(pngdata,"plot.png");
}


//linear
void linreg::calcmata()
{
    //cout<<"mata";
    mata[0][0]=num;
    mata[0][1]=mata[1][0]=sum(x,num,1);
    mata[1][1]=sum(x,num,2);
}

void linreg::calcmatb()
{
    //cout<<"matb";
    matb[0]=sum(y,num,1);
    matb[1]=0;
    for(int i=0;i<num;i++)
    {
        matb[1]+=x[i]*y[i];
    }
}

void linreg::calcparam()
{
    //cout<<"calcparam";
   float ainv[5][5];
   float det=mata[0][0]*mata[1][1]-mata[0][1]*mata[1][0];
   ostringstream str11;


   ainv[0][0]=mata[1][1]/det;
   ainv[1][1]=mata[0][0]/det;
   ainv[0][1]=(-1*mata[0][1])/det;
   ainv[1][0]=(-1*mata[1][0])/det;

   b[0]=ainv[0][0]*matb[0]+ainv[0][1]*matb[1];
   b[1]=ainv[1][0]*matb[0]+ainv[1][1]*matb[1];

   str11<<fixed;
    if(b[1]==1)
        str11<<"x ";
    else if(b[1]==-1)
        str11<<"-x ";
    else
        str11<<setprecision(2)<<b[1]<<"x ";
    if(b[0]>0)
        str11<<"+ ";
    if(b[0]!=0)
        str11<<setprecision(2)<<b[0];

    exp.append(str11.str());
    cout<<exp<<endl;
        //cout<<ainv[0][0]<<" "<<ainv[0][1]<<" "<<ainv[1][0]<<" "<<ainv[1][1]<<"\n";
   //cout<<"\n"<<b[1]<<"x+ "<<b[0]<<"\n";
}

void linreg::calcpredy()
{
    int cnt=0;
    for(double i=minarr(x,num)-5;i<=maxarr(x,num)+5;i+=0.1)
    {
        ypred[cnt]=(b[1]*i)+b[0];
       // cout<<ypred[cnt]<<"\n";
        cnt++;
    }
}

//quadratic
void quadreg::calcmata()
{
    mata[0][0]=num;
    mata[1][0]=mata[0][1]=sum(x,num,1);
    mata[2][0]=mata[0][2]=mata[1][1]=sum(x,num,2);
    mata[1][2]=mata[2][1]=sum(x,num,3);
    mata[2][2]=sum(x,num,4);

    /*for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
            cout<<mata[i][j]<<" ";
        cout<<"\n";
    }*/
}

void quadreg::calcmatb()
{
    matb[0]=sum(y,num,1);
    matb[1]=0;
    for(int i=0;i<num;i++)
    {
        matb[1]+=x[i]*y[i];
    }
    matb[2]=0;
    for(int i=0;i<num;i++)
    {
        matb[2]+=x[i]*x[i]*y[i];
    }

    //cout<<"\n"<<matb[0]<<" "<<matb[1]<<" "<<matb[2]<<"\n";
}

void quadreg::calcparam()
{
    float det=0,ainv[5][5];
    for(int i=0;i<3;i++)
    {
        det=det+(mata[0][i]*(mata[1][(i+1)%3]*mata[2][(i+2)%3] - mata[1][(i+2)%3]*mata[2][(i+1)%3]));
    }

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            ainv[i][j]=((mata[(i+1)%3][(j+1)%3] * mata[(i+2)%3][(j+2)%3]) - (mata[(i+1)%3][(j+2)%3]* mata[(i+2)%3][(j+1)%3]))/det;
        }
    }

    b[0]=ainv[0][0]*matb[0]+ainv[0][1]*matb[1]+ainv[0][2]*matb[2];
    b[1]=ainv[1][0]*matb[0]+ainv[1][1]*matb[1]+ainv[1][2]*matb[2];
    b[2]=ainv[2][0]*matb[0]+ainv[2][1]*matb[1]+ainv[2][2]*matb[2];

    /*for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
            cout<<ainv[i][j]<<" ";
        cout<<"\n";
    }*/
    ostringstream str11;
    str11<<fixed;
    for(int i=2;i>0;i--)
    {
        b[i]=round(b[i]*100)/100.00;
        if(i!=2 && b[i]>0)
        {
            str11<<"+ ";
        }
        if(b[i]==1)
            str11<<"x^"<<i<<" ";
        else if(b[i]==-1)
            str11<<"-x^"<<i<<" ";
        else
        {
            if(b[i]!=0)
                str11<<setprecision(2)<<b[i]<<"x^"<<i<<" ";
        }
    }
    if(b[0]!=0)
    {
        if(b[0]>0)
            str11<<"+ ";
        str11<<setprecision(2)<<b[0];
    }
    exp.append(str11.str());
    cout<<exp<<endl;

    //cout<<b[2]<<"x^2+ "<<b[1]<<"x+ "<<b[0];
}

void quadreg::calcpredy()
{
    int cnt=0;
    for(double i=minarr(x,num)-5;i<=maxarr(x,num)+5;i+=0.1)
    {
        ypred[cnt]=(b[2]*i*i)+(b[1]*i)+b[0];
        //cout<<ypred[cnt];
        cnt++;
    }
}

//cubic
void cubreg::calcmata()
{
    mata[0][0]=num;
    mata[1][0]=mata[0][1]=sum(x,num,1);
    mata[2][0]=mata[1][1]=mata[0][2]=sum(x,num,2);
    mata[3][0]=mata[2][1]=mata[1][2]=mata[0][3]=sum(x,num,3);
    mata[3][1]=mata[2][2]=mata[1][3]=sum(x,num,4);
    mata[3][2]=mata[2][3]=sum(x,num,5);
    mata[3][3]=sum(x,num,6);
}

void cubreg::calcmatb()
{
    matb[0]=sum(y,num,1);
    matb[1]=0;
    for(int i=0;i<num;i++)
    {
        matb[1]+=x[i]*y[i];
    }
    matb[2]=0;
    for(int i=0;i<num;i++)
    {
        matb[2]+=x[i]*x[i]*y[i];
    }
    matb[3]=0;
    for(int i=0;i<num;i++)
    {
        matb[3]+=pow(x[i],3)*y[i];
    }
}

void cubreg::calcparam()
{
    float ainv[5][5];
    inverse(mata,ainv);

    for(int i=0;i<4;i++)
    {
        b[i]=0;
        for(int j=0;j<4;j++)
        {
            b[i]+=ainv[i][j]*matb[j];
        }
    }

    ostringstream str11;
    str11<<fixed;
    for(int i=3;i>0;i--)
    {
        b[i]=round(b[i]*100)/100.00;
        if(i!=3 && b[i]>0)
        {
            str11<<"+ ";
        }
        if(b[i]==1)
            str11<<"x^"<<i<<" ";
        else if(b[i]==-1)
            str11<<"-x^"<<i<<" ";
        else
        {
            if(b[i]!=0)
                str11<<setprecision(2)<<b[i]<<"x^"<<i<<" ";
        }
    }
    if(b[0]!=0)
    {
        if(b[0]>0)
            str11<<"+ ";
        str11<<setprecision(2)<<b[0];
    }
    exp.append(str11.str());
    cout<<exp<<"\n";

    //cout<<b[3]<<" "<<b[2]<<" "<<b[1]<<" "<<b[0]<<"\n";
}

void cubreg::calcpredy()
{
    int cnt=0;
    for(double i=minarr(x,num)-5;i<=maxarr(x,num)+5;i+=0.1)
    {
        ypred[cnt]=(b[3]*i*i*i)+(b[2]*i*i)+(b[1]*i)+b[0];
        cnt++;
    }
}

