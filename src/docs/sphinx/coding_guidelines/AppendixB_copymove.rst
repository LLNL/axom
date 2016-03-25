************************************
Appendix B Copy and Move Operations
************************************

This appendix discusses proper implementation of copy operations (copy 
constructor and copy-assignment operator) and move operations (move 
constructor and move-assignment operator). C++ copy semantics have been 
part of the language for many years and are often required when using STL 
containers, supporting value semantics, etc. 

Move semantics are a feature introduced in C++11. Before mover semantics 
were introduced, copying was the only way to transfer state (i.e., the 
collective set of non-static data member values) from one object
to another. Copying causes the destination object to have the same state as 
the source and does not modify the source. Moving causes the destination
object to acquire the state of the source and leaves the source object in 
an **unspecified** state. 

Move operations are not required since copy operations can be used to 
accomplish the same result. Rather, move semantics are an important 
optimization feature of C++11, whereby relatively expensive copy operations 
may be replaced by cheaper *move* operations. Also, the C++11 standard 
requires compilers to use move operations instead of copy operations when 
certain conditions are fulfilled.

Often, compiler-generated methods (described in Appendix A) are insufficient
for classes that manage resources. In the next couple of sections, we describe 
proper implementation of copy and move operations using the following simple
array class as a concrete example:

.. code-block:: cpp

   template <typename T>
   class MyArray
   {
   public:

      // Standard constructors and destructor omitted for brevity.
      // However, it is important to emphasize that compiler-generated
      // versions of these methods are insufficient since the class
      // manages resources, namely its array of data.

      //
      // Copy constructor.
      // 
      MyArray(const MyArray<T>& other);

      //
      // Copy-assignment operator.
      // 
      MyArray<T>& operator=(const MyArray<T>& rhs);

      //
      // Move constructor.
      // 
      MyArray(MyArray<T>&& other);

      //
      // Move-assignment operator.
      // 
      MyArray<T>& operator=(MyArray<T>&& rhs);

      //
      // Swap method to use the copy-and-swap idiom.
      // 
      void swap(MyArray& other);

      // other functions...

   private:
      T* m_data;
      size_t m_size;
   };


========================================================
B.1 Copy operations 
========================================================

This section shows how to implement a copy constructor and a 
copy-assignment operator for the MyArray class in a traditional
manner found in many C++ codes. Alternative implementations are
possible and should be pursued if it is necessary to make these
methods *exception-safe*, for example. 

The copy constructor is implemented as follows:

.. code-block:: cpp

   template <typename T>
   MyArray<T>::MyArray(const MyArray<T>& other)
   : m_data(new T[other.m_size]),
     m_size(other.m_size)
   {
      std::copy(other.m_data, other.m_data + m_size, m_data);
   }

Note that we use the standard library function to copy data from the
input object to *this* object. A raw for-loop could be used instead and
is an acceptable alternative.

The copy-assignement operator is implemented as follows:

.. code-block:: cpp

   template <typename T>
   MyArray<T>& MyArray<T>::operator=(const MyArray<T>& rhs)
   {
      if (this != &rhs) {
         MyArray<T> temp(rhs);
         swap(rhs);
      }
      return *this;
   }

First, we check for assignment to self to avoid temporary object creation 
and copy operations when unnecessary. Then, we use the *copy-and-swap* idiom 
by calling the copy constructor to create a temporary object and then use the 
swap method provided by the class to copy the class data members from the 
temporary object. Finally, we return a reference to "\*this". 

The swap method for the class is implemented as:

.. code-block:: cpp

   template <typename T>
   void MyArray<T>::swap(MyArray<T>& rhs)
   {
      std::swap(m_data, other.m_data);
      std::swap(m_size, other.m_size);
   }


========================================================
B.2 Move operations 
========================================================

This section shows how to implement a move constructor and a 
move-assignment operator for the MyArray class.

Recall that earlier we said that a move operation causes the destination object 
to acquire the state of the source and leaves the source object in an 
**unspecified** state. Thus, one should always assume that the source object
no longer owns any resources and that its state is similar to an empty
object. In other words, a move operation does not allocate new resources,
as a copy operation does. Instead it steals resources from one object and
gives them to another.

Thus, the move constructor is implemented as follows:

.. code-block:: cpp

   template <typename T>
   MyArray<T>::MyArray(MyArray<T>&& other)
   {
      m_data = other.m_data;
      m_size = other.m_size;

      other.m_data = ATK_NULLPTR;
      other.m_size = 0;
   }

Note that the argument is a (non-const) C++11 rvalue reference variable, or 
*universal reference*. The constructor simply moves the state from the 
argument object to "this" object. Since it does not allocate memory or copy
any memory buffers, the move constructor is potentially much faster than the
copy constructor.

The move-assignement operator is implemented as follows:

.. code-block:: cpp

   template <typename T>
   MyArray<T>& MyArray<T>::operator=(MyArray<T>&& rhs)
   {      
      if (this != &rhs) {

         delete [ ] m_data;
         m_size = 0;

         m_data = other.m_data;
         m_size = other.m_size;

         other.m_data = ATK_NULLPTR;
         other.m_size = 0;

      }
      return *this;
   }

First, we check for assignment to self to avoid potentially unwanted data
destruction and invalidation of object state. Then, we release any resources
"this" object owns, take ownership of the other object's resources, and set
the other object state to be "empty".  Finally, we return a reference to 
"\*this".


========================================================
B.3 Overloaded move and copy operations
========================================================

C++11 overload resolution methods were modified to support *rvalue references*
or *universal references* (briefly described in Appendix B). Technically,
an rvalue is an unnamed value that exists only during the evaluation of an
expression. 

What does this mean?

For example, the standard library vector class method "vector::push_back()"
now has two overloaded versions. One takes a parameter of type "const T&" 
for *lvalue* arguments as before. The new one takes an *rvalue* reference
of type "T&&". Which method is called depends on usage. 

For example, recall the "MyArray" class defined above. The following code 
example will add two MyArray objects to a standard vector:

.. code-block:: cpp
  
    #include <vector>
    
    int main()
    {
       std::vector< MyArray<int> > vec;
       vec.push_back( MyArray<int>(10) );
       vec.push_back( MyArray<int>(20) );
    } 

Here, both calls to "push_back" resolve to "push_back(T&&)" because their
arguments are *rvalues*. These operations *move* the arguments' resources into 
the objects in the vector using the MyArray move constructor. Prior to C++11,
the same code would generate copies of the argument objects using the copy
constructor.

The following code example will resolve to the "push_back(const T&)" version
and perform copy operations since the argument is an *lvalue*:

.. code-block:: cpp
 
    #include <vector>
   
    int main()
    {
       std::vector< MyArray<int> > vec;
       MyArray<int> ma(10);
       vec.push_back( ma );
    }

This "push_back" operation uses the MyArray copy constructor.

Incidentally, the selection of "push_back(T&&)" can be forced by casting the
lvalue to an rvalue reference::

    vec.push_back( static_cast< MyArray<int>&& >(ma) );

Alternatively, one could use the new standard function "std::move()" for
the same purpose; i.e.::

    vec.push_back( std::move(ma) );

While move operations are often the best choice because they eliminate 
unnecessary copies, the copy constructor and copy-assignement operator
are necessary to retain pure copy semantics -- when we want an argument 
to "std::vector::push_back()" to retain its state, for example.
