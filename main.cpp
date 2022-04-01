#include "calc_classdefn.h"

#include "pbPlots.hpp"
#include "supportLib.hpp"
#include<iostream>
#include<iomanip>
#include<string>
#include<stdlib.h>
#include<math.h>
#define PI 3.14


using namespace std;

int main()
{
    char c;
    int ch,select;
    expeval e;
    comp c1,c2;
    quadroot q;
    reg *r;

while(1)
{
    system("CLS");
    cout<<"\n\n\t\t\t*************************************\n\n\t\t\tSCIENTIFIC CALCULATOR WITH REGRESSION\n";
    cout<<"\n\t\t\t*************************************\n\n";

    cout<<"\t\t\t1. Calculate (Basic,trig,log)\n\n\t\t\t2. Complex Arithmetic\n\n\t\t\t3. Quadratic root finder\n\n\t\t\t4. Regression & Graphing\n\n";
    cout<<"\n\t\t\tEnter choice(Press 0 to exit): ";
    cin>>select;
    system("CLS");
    switch(select)
    {
    case 0:
        exit(0);

    case 1:
        cout<<"\n\n\t\t\t**********************************\n\n\t\t\tSCIENTIFIC CALCULATOR INTERPRETER\n";
        cout<<"\n\t\t\t*********************************\n\n";
        e.dispinst();
        e.run();
        break;

    case 2:
            do
            {
                system("CLS");
                cout<<"\n\n\t\t\t********************\n\n\t\t\t COMPLEX CALCULATOR\n";
                cout<<"\n\t\t\t********************\n\n";
                cout<<"\nFor first number\n";
                cin>>c1;
                cout<<"\nFor second number\n";
                cin>>c2;
                ch=compmenu(c1,c2);

                cout<<"\n\nDo you wish to continue(Y/N) ? ";
                cin>>c;
            }while(c=='y'||c=='Y');
            break;

    case 3:

        do
        {
            system("CLS");
            cout<<"\n\n\t\t\t**********************\n\n\t\t\tQUADRATIC ROOT FINDER\n";
            cout<<"\n\t\t\t**********************\n\n";
            cout<<"Enter polynomial:\n\n";
            cin>>q;
            q.compute();
            q.disppoly();
            cout<<q;

            cout<<"\n\nDo you wish to continue (Y/N) ? ";
            cin>>c;
        }while(c=='Y'||c=='y');
        break;

    case 4:
        do
        {
        system("CLS");
        cout<<"\n\n\t\t\t*******************************\n\n\t\t\tREGRESSION CALCULATOR & PLOTTER\n";
        cout<<"\n\t\t\t*******************************\n\n";
        cout<<"\n\nDo you wish to perform:\n\t\t1. Linear regression\n\t\t2. Quadratic regression\n\t\t3. Cubic regression\n\n\t\tEnter choice: ";
        cin>>ch;

        switch(ch)
        {
        case 1:
            r=new linreg;
            break;
        case 2:
            r=new quadreg;
            break;
        case 3:
            r=new cubreg;
            break;
        }

        r->input();
        r->disp();
        r->calcmata();
        r->calcmatb();
        cout<<"\nThe required equation is:\n\t\ty = ";
        r->calcparam();
        r->calcpredy();
        r->plot();
        cout<<"\n\nThe plot has been saved in the current directory.";
        cout<<"\n";

        cout<<"\n\nDo you wish to continue (Y/N) ? ";
        cin>>c;
        }while(c=='y'||c=='Y');

        break;
    }
}
    return 0;
}
