#ifndef VTENSE_HPP
#define VTENSE_HPP
#include "point.hpp"

namespace geom {
  template <class vtype,class CRefInner,int d>
    struct cref {
       cref(CRefInner  cr_):cr(cr_){}
       
       operator vtype() const { 
          return vtype(cr);
       }

       cref& operator=(const vtype &v)
       {
           cr = geom::eigen_vec(v);
           return *this;
       }
       
       auto operator[](size_t sz) const -> typename std::remove_reference<decltype(vtype()[sz])>::type { 
           return cr[sz];
       }
       auto operator-(const vtype& o) ->  decltype( vtype() - vtype())
       {
           return vtype(*this) - o;
       }
       
       template <class Vec>
       cref & operator +=(const Vec& o) 
       {
           cr += geom::eigen_vec(o);
           return *this;
       }
       template <class Vec>
       cref & operator -=(const Vec& o) 
       {
           cr -= geom::eigen_vec(o);
           return *this;
       }
    private:
       CRefInner cr;
    };
    template <class vtype,class CRefInner>
    struct cref<vtype, CRefInner,1> {
      
       cref(CRefInner  cr_):cr(cr_){}

       operator vtype() const { 
          return cr[0];
       }
       cref& operator=(const vtype &v)
       {
           cr[0] = v;
           return *this;
       }
    private:
       CRefInner cr;
    };

template <class vtype>
struct vector_of_type
{
    enum{dim = geom::point_dim < vtype >::dimension};

    typedef Eigen::Matrix<double,Eigen::Dynamic,dim> EMT;
    typedef vtype value_type;
    typedef Eigen::Block<EMT,1,dim,false> value_type_inner;
    typedef const Eigen::Block<const EMT,1,dim,false>& const_reference_inner;
    typedef Eigen::Block<EMT,1,dim,false>& reference_inner;

  
    typedef cref<vtype,const_reference_inner,dim> const_reference;
    typedef cref<vtype,reference_inner,dim> reference;

    vector_of_type():sz(0){}
    vector_of_type(int sz_):sz(sz_){
        reserve(sz);
    }

    vector_of_type(vector_of_type && o)
    {
        cont.swap(o.cont);
    }

    vector_of_type(const vector_of_type & o):cont(o.cont)
    {
    }
  
    template <class Iter>
    vector_of_type(Iter f, Iter l)
    {
        reserve(std::distance(f,l));
        std::copy(f,l,std::back_inserter(*this));
    }
    template <class Iter>
    void assign(Iter f, Iter l)
    {
        std::copy(f,l,std::back_inserter(*this));
    }
    const_reference operator[](int i) const {
        return cont.row(i);
    }
   
    reference operator[](int i) {
        return cont.row(i);
    }

    void reserve(size_t sz) {
        cont.resize(Eigen::NoChange_t(),sz);
    }

    void push_back(const vtype& p) {
        cont.row(sz) = geom::eigen_vec(p);
        ++sz;
    }

    const_reference back() const  { 
        return cont.row(sz-1);
    }
    reference back() { 
        return cont.row(sz-1);
    }

    reference front() { 
        return cont.row(0);
    }
    const_reference front() const { 
        return cont.row(0);
    }
    
    
    struct const_row_iterator : std::iterator<std::random_access_iterator_tag,value_type,ptrdiff_t,value_type*,const_reference> {
        
        typedef vector_of_type container_type;

        const_row_iterator (const const_row_iterator& o):cont(o.cont),row(o.row) {
        }
        const_row_iterator (const vector_of_type& c,int row_) : cont(c),row(row_){ 
        }
        const_row_iterator& operator=(const const_row_iterator&o ) 
        {
            row = o.row;
            return *this;
        }
        const_reference operator*() const { return cont[row]; }
        const_reference operator[](size_t sz) const { return cont[row+sz]; }
        const_row_iterator& operator++() { ++row; return *this;}
        const_row_iterator operator++(int) { ++row; return const_row_iterator(cont,row-1);}
        
        const_row_iterator& operator--() { --row; return *this;}
        const_row_iterator operator--(int) { --row; return const_row_iterator(cont,row+1);}
        
        const_row_iterator operator-(difference_type sz) const
        {
          return const_row_iterator(cont, row - sz);
        }
        ptrdiff_t operator-(const_row_iterator o) const
        {
          assert(&cont==&o.cont);
          return row - o.row;
        }
        const_row_iterator &operator+=(difference_type sz) { 
            row+=sz;
            return *this;
        }

        const_row_iterator operator+(difference_type sz) const
        {
          return const_row_iterator(cont, row + sz);
        }
        bool operator==(const const_row_iterator& o) const { return &o.cont==&cont && row == o.row; }
        bool operator!=(const const_row_iterator& o) const { return !(o == *this); }
        bool operator<(const const_row_iterator& o) const { return &o.cont==&cont && row < o.row; }

        int row;
        
        const vector_of_type& cont;
    };

    struct const_rev_row_iterator : std::iterator<std::random_access_iterator_tag,value_type,ptrdiff_t,value_type*,const_reference> {
        
         typedef vector_of_type container_type;
        const_rev_row_iterator (const const_rev_row_iterator& o):cont(o.cont),row(o.row){
        }
        const_rev_row_iterator (const vector_of_type& c,int row_):cont(c),row(row_){ 
        }
        const_rev_row_iterator& operator=(const const_rev_row_iterator&o ) 
        {
            row = o.row;
            return *this;
        }
        const_reference operator*() const { return  (cont[row]); }
        const_reference operator[](size_t sz) const { return (cont[row]); }
        const_rev_row_iterator operator++() { --row; return *this;}
        const_rev_row_iterator operator++(int) { --row; return const_rev_row_iterator(cont,row+1);}

        const_rev_row_iterator operator--() { ++row; return *this;}
        const_rev_row_iterator operator--(int) { ++row; return const_rev_row_iterator(cont,row-1);}


        const_rev_row_iterator operator-(difference_type sz) const
        {
          return const_rev_row_iterator(cont, row - sz);
        }
        const_rev_row_iterator operator+(difference_type sz) const
        {
          return const_rev_row_iterator(cont, row + sz);
        }
        ptrdiff_t operator-(const_rev_row_iterator o) const
        {
          assert(&cont==&o.cont);
          return o.row - row;
        }
        const_rev_row_iterator &operator+=(difference_type sz) { 
            row-=sz;
            return *this;
        }
         bool operator==(const const_rev_row_iterator& o) const { return &o.cont==&cont && row == o.row; }
         bool operator!=(const const_rev_row_iterator& o) const { return !(o==*this); }
         bool operator<(const const_rev_row_iterator& o) const { return &o.cont==&cont && row > o.row; }
         int row;
        const vector_of_type& cont;
    };

    struct row_iterator : std::iterator<std::random_access_iterator_tag,value_type,ptrdiff_t,value_type*,reference> {
        
        typedef vector_of_type container_type;
        row_iterator (row_iterator& o):cont(o.cont),row(o.row) {
        }

        row_iterator (vector_of_type& c,int row_):cont(c),row(row_) { 
        }

        row_iterator& operator=(const row_iterator&o ) 
        {
            row = o.row;
            return *this;
        }

        row_iterator operator-(difference_type sz) const
        {
          return row_iterator(cont, row - sz);
        }
        row_iterator operator+(difference_type sz) const
        {
          return row_iterator(cont, row + sz);
        }
        ptrdiff_t operator-(row_iterator o) const
        {
          assert(&cont==&o.cont);
          return row - o.row;
        }
         row_iterator &operator+=(difference_type sz) { 
            row+=sz;
            return *this;
        }
       
        reference operator*()             { return cont[row]; }
        reference operator[](size_t sz)   { return cont[row+sz]; }

        row_iterator operator++()    { ++row; return *this;}
        row_iterator operator++(int) { ++row; return row_iterator(cont, row-1);}
        
        row_iterator& operator--() { --row; return *this;}
        row_iterator operator--(int) { --row; return row_iterator(cont,row+1);}
 
        bool operator==(const row_iterator& o) const{ return &o.cont==&cont && row == o.row; }
        bool operator!=(const row_iterator& o) const { return !(o==*this); }
        bool operator<(const row_iterator& o) const { return &o.cont==&cont && row < o.row; }
        int row;
        vector_of_type& cont;
    };
    
     struct rev_row_iterator : std::iterator<std::random_access_iterator_tag,value_type,ptrdiff_t,value_type*,reference> {
                
        typedef vector_of_type container_type;
        rev_row_iterator(const row_iterator& o)
            :cont(o.cont),row(o.row) {
        }

        rev_row_iterator (vector_of_type& c,int row_):cont(c),row(row_) { 
        }
        rev_row_iterator& operator=(const rev_row_iterator&o ) 
        {
            row = o.row;
            return *this;
        }

        reference operator*()             { return cont[row]; }
        reference operator[](size_t sz)   { return cont[row+sz]; }

        rev_row_iterator operator++()    { return --row;return *this;}
        rev_row_iterator operator++(int) { return --row; return rev_row_iterator(cont, row+1);}
        
        rev_row_iterator operator--() { ++row; return *this;}
        rev_row_iterator operator--(int) { ++row; return rev_row_iterator(cont,row-1);}


        rev_row_iterator operator-(difference_type sz) const
        {
          return rev_row_iterator(cont, row - sz);
        }
        rev_row_iterator operator+(difference_type sz) const
        {
          return rev_row_iterator(cont, row + sz);
        }
        ptrdiff_t operator-(rev_row_iterator o) const
        {
          assert(&cont==&o.cont);
          return o.row - row;
        }

         rev_row_iterator &operator+=(difference_type sz) { 
            row-=sz;
            return *this;
        }

        bool operator==(const rev_row_iterator& o) const{ return &o.cont==&cont && row == o.row; }
        bool operator!=(const rev_row_iterator& o) const { return !(o==*this); }
        bool operator<(const rev_row_iterator& o) const { return &o.cont==&cont && row > o.row; }
        int row;
        vector_of_type& cont;
    };
    
    const_row_iterator  cbegin() const { return const_row_iterator (*this,0);         }
    const_row_iterator  cend()   const { return const_row_iterator (*this,sz); }

    const_row_iterator  begin() const { return const_row_iterator (*this,0);         }
    const_row_iterator  end()   const { return const_row_iterator (*this,sz); }

    const_rev_row_iterator crbegin() const { return const_rev_row_iterator(*this, size() -1); }
    const_rev_row_iterator crend()   const { return const_rev_row_iterator(*this, -1); }


    row_iterator  begin() { return row_iterator (*this,0);         }
    row_iterator  end()   { return row_iterator (*this,sz); }

    rev_row_iterator rbegin() { return rev_row_iterator(*this, size() -1); }
    rev_row_iterator rend()   { return rev_row_iterator(*this, -1); }

    void swap(vector_of_type & o)
    {
        cont.swap(o.cont);
    }
    size_t size() const { return sz;}

    EMT& get() { return cont; }
    const EMT& cget() const { return cont; }
    EMT cont;
    size_t sz;
};


/* template <class vtype>
struct matrix_of_type
{
    enum{dim = point_dim < vtype >::dimension};

    typedef Eigen::MatrixXd EMT;
    typedef vtype value_type;

    matrix_of_type(){}
    matrix_of_type(int szm, int szn)
    {
        for(int i = 0;i < dim; ++i)
            cont.resize(szm, szn);
    }
    matrix_of_type(matrix_of_type && o)
    {
        swap(o);
    }
    matrix_of_type(const matrix_of_type & o)
    {
        for(int i = 0;i < dim; ++i)
            cont[i] = o.cont[i];
    }

    vtype operator()(int i, int j) const {
        return vtype (cont(i *sizem + j));
    }

    void reserve(int szm, int szn) {
        for(int i = 0;i < dim; ++i)
            cont.resize(szm, szn);
    }

    void push_back(pt_t < dim >& p) {
        cont[size] = eigen_vec(p);
        ++size;
    }

    void swap(matrix_of_type & o)
    {
        for(int i = 0;i < dim; ++i)
            cont[i].swap(o.cont[i]);
    }

    EMT& get() { return cont; }
    const EMT& cget() const { return cont; }
    EMT cont[dim];
}; */ 
}
    
     namespace std {

          template <class vtype, class T,int d>
     struct  remove_reference<  geom::cref<vtype,T,d> >
     {
         typedef vtype type;
     };
template<class vtype>
struct _Is_checked_helper< decltype(geom::vector_of_type<vtype>().begin()) >
    : public true_type
{   // mark row_iterator as checked
};


template<class vtype>
struct _Is_checked_helper< decltype(geom::vector_of_type<vtype>().cbegin()) >
    : public true_type
{   // mark const_row_iterator as checked
};



template<class vtype>
struct _Is_checked_helper< decltype(geom::vector_of_type<vtype>().rbegin()) >
    : public true_type
{   // mark rev_row_iterator as checked
};


template<class vtype>
struct _Is_checked_helper< decltype(geom::vector_of_type<vtype>().crbegin()) >
    : public true_type
{   // mark rev_const_row_iterator as checked
};

}


#endif // VTENSE_HPP
