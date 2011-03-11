// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MULTIVECTOR_HH
#define DUNE_MULTIVECTOR_HH

#include <dune/common/iteratorfacades.hh>
#include <dune/common/typetraits.hh>
#include <vector>

#include <iostream>
#include <cassert>

namespace Dune {

  /**
     proxy object to give access
   */
  template<typename A, typename B, typename C, typename D>
  struct MultiDataProxy;

  template<typename A, typename B, typename C, typename D>
  struct MultiDataProxy
  {
    typedef MultiDataProxy<typename remove_const<A>::type,
        typename remove_const<B>::type,
        typename remove_const<C>::type,
        typename remove_const<D>::type> MutableProxy;
    typedef MultiDataProxy<const typename remove_const<A>::type,
        const typename remove_const<B>::type,
        const typename remove_const<C>::type,
        const typename remove_const<D>::type> ConstProxy;

    A & _a;
    B & _b;
    C & _c;
    D & _d;

    int pos;
    std::string name;
    MultiDataProxy(A& a, B& b, C& c, D& d, size_t pos, std::string _n) :
      _a(a), _b(b), _c(c), _d(d), pos(pos), name(_n) {}
    MultiDataProxy(ConstProxy & other) :
      _a(other._a), _b(other._b), _c(other._c), _d(other._d), pos(other.pos), name(other.name) {}
    MultiDataProxy(MutableProxy & other) :
      _a(other._a), _b(other._b), _c(other._c), _d(other._d), pos(other.pos), name(other.name) {}
    // compare
    bool operator == (const MultiDataProxy & other) const { return a() == other.a(); }
    bool operator != (const MultiDataProxy & other) const { return a() != other.a(); }
    bool operator  < (const MultiDataProxy & other) const { return a()  < other.a(); }
    bool operator  > (const MultiDataProxy & other) const { return a()  > other.a(); }
    // assign
    MultiDataProxy& operator = (const ConstProxy & other) {
      assign(other);
      return *this;
    }
    MultiDataProxy& operator = (const MutableProxy & other) {
      assign(other);
      return *this;
    }

    // access
    typename A::reference  a () { return _a[pos]; }
    typename B::reference  b () { return _b[pos]; }
    typename C::reference  c () { return _c[pos]; }
    typename D::reference  d () { return _d[pos]; }
    typename A::value_type a () const { return _a[pos]; }
    typename B::value_type b () const { return _b[pos]; }
    typename C::value_type c () const { return _c[pos]; }
    typename D::value_type d () const { return _d[pos]; }
  private:
    template<typename P>
    void assign(const P & other) {
#ifdef DEBUG_MULTIVEC
      std::cerr << "Assign " << name << "," << pos
                << "\n   from " << other.name << "," << other.pos << std::endl;
#endif
      a() = other.a();
      b() = other.b();
      c() = other.c();
      d() = other.d();
    }
  };

  template<typename A, typename B, typename C, typename D>
  std::ostream& operator<< (std::ostream & s, const MultiDataProxy<A,B,C,D> & i)
  {
    return s << "("
           << i.a() << ", "
           << i.b() << ", "
           << i.c() << ", "
           << i.d() << ")";
  }

  template<typename A, typename B, typename C, typename D>
  class MultiVectorIterator :
    public Dune::BidirectionalIteratorFacade< MultiVectorIterator<A,B,C,D>,
        MultiDataProxy<A,B,C,D> >
  {
    // friend class MultiVectorIterator<typename remove_const<C>::type, typename remove_const<T>::type >;
    // friend class TestIterator<const typename remove_const<C>::type, const typename remove_const<T>::type >;
  public:
    mutable MultiDataProxy<A,B,C,D> data;

    // constructors
    MultiVectorIterator(A& a, B& b, C& c, D& d, size_t n, std::string i) :
      data(a,b,c,d, n, i) {}
    MultiVectorIterator(const MultiVectorIterator & other) :
      data(other.data)
    {
#ifdef DEBUG_MULTIVEC
      std::cerr << "Copy Iterator " << data.name << "," << data.pos << std::endl;
#endif

    }

    size_t pos() const { return data.pos; }

    MultiVectorIterator operator = (const MultiVectorIterator & other)
    {
#ifdef DEBUG_MULTIVEC
      // std::cerr << "Assign Iterator " << data.name << "," << data.pos
      //           << "\n   from " << other.data.name << "," << other.data.pos << std::endl;
#endif
      assert(other.data._a == data._a
             && other.data._b == data._b
             && other.data._c == data._c
             && other.data._d == data._d);

      data.pos = other.data.pos;
      return *this;
    }

    // operators
    bool equals (const MultiVectorIterator & other) const
    {
#ifdef DEBUG_MULTIVEC
      // std::cerr << "Compare " << data.name << "," << data.pos
      //           << " with " << other.data.name << "," << other.data.pos << "\n";
#endif
      assert(other.data._a == data._a
             && other.data._b == data._b
             && other.data._c == data._c
             && other.data._d == data._d);

      return other.data.pos == data.pos;
      // && other.data._a == data._a
      // && other.data._b == data._b
      // && other.data._c == data._c
      // && other.data._d == data._d;
    }
    // in-/decrement
    void increment()
    {
#ifdef DEBUG_MULTIVEC
      // std::cerr << "Increment " << data.name << "," << data.pos << std::endl;
#endif
      data.pos++;
    }
    void decrement()
    {
#ifdef DEBUG_MULTIVEC
      std::cerr << "Decrement " << data.name << "," << data.pos << std::endl;
#endif
      data.pos--;
    }
    // dereference
    MultiDataProxy<A,B,C,D>& dereference() const
    {
#ifdef DEBUG_MULTIVEC
      std::cerr << "dereference " << data.name << "," << data.pos << std::endl;
#endif
      return data;
    }
  };

  template<typename A, typename B, typename C, typename D>
  class MultiVector
  {
    typedef MultiVector<A,B,C,D> vector_type;

    A & _a;
    B & _b;
    C & _c;
    D & _d;

    std::string name;

  public:
    //! container interface typedefs
    //! \{

    /** \brief Type of the values stored by the container */
    typedef MultiDataProxy<A,B,C,D> value_type;

    /** \brief Type of the const values stored by the container */
    typedef MultiDataProxy<A,B,C,D> const_value_type;

    /** \brief Reference to a small block of bits */
    typedef value_type reference;

    /** \brief Const reference to a small block of bits */
    typedef const_value_type const_reference;

    /** \brief Pointer to a small block of bits */
    typedef reference* pointer;

    /** \brief Const pointer to a small block of bits */
    typedef const_reference* const_pointer;

    /** \brief size type */
    typedef size_t size_type;

    //! \}

    //! iterators
    //! \{
    typedef MultiVectorIterator<A,B,C,D> iterator;
    typedef MultiVectorIterator<const A, const B, const C, const D> const_iterator;
    //! \}

    MultiVector(A & a, B & b, C & c, D & d, std::string n = "?") :
      _a(a), _b(b), _c(c), _d(d)
    {
      assertSize();
      name = n;
    }

    iterator begin()
    {
      return iterator(_a,_b,_c,_d, 0, name);
    }

    iterator end()
    {
      assertSize();
      return iterator(_a,_b,_c,_d, _a.size(), name);
    }

    void erase(iterator front, iterator back)
    {
      assertSize();
      _a.erase( _a.begin()+front.pos(), _a.begin()+back.pos());
      _b.erase( _b.begin()+front.pos(), _b.begin()+back.pos());
      _c.erase( _c.begin()+front.pos(), _c.begin()+back.pos());
      _d.erase( _d.begin()+front.pos(), _d.begin()+back.pos());
    }
    size_t size() const {
      assertSize();
      return _a.size();
    }
  private:
    void assertSize() const {
      assert(_a.size() == _b.size() &&
             _a.size() == _c.size() &&
             _a.size() == _d.size() );
    }
  };

} // end namespace Dune

#endif // DUNE_MULTIVECTOR_HH
