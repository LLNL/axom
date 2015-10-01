
/* Compute factorial of n */
int  fact(int n);

int fact(int n)
{
  if (n > 1)
    return n*fact(n-1);
  else
    return 1;
}

typedef struct FooStruct
{
   int a;
   double b;
} Foo;

class Darray
{
 public:
   Darray(int size=10){n=size; d = new double[n];}
   ~Darray(){delete [] d;}
   int n;
   double * d;
   double get(int i) { return d[i]; }
   void set(int i, double x) { d[i] = x; }
   int size() {return n;}
   
};
