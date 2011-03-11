// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MULTIVECTOR_HH
#define DUNE_MULTIVECTOR_HH

#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/typetraits.hh>
#include <vector>

#include <iostream>
#include <cassert>

namespace Dune {

#ifndef DOXYGEN
  // TMPs for the MultiVector tuples stuff
  namespace {
    struct MultiVectorAssertSize
    {
      size_t sz_;
      MultiVectorAssertSize(size_t sz) : sz_(sz) {}
      template <class E>
      void visit(E& elem) { assert(elem.size() == sz_); }
    };
    struct MultiVectorEraser
    {
      size_t front, back;
      MultiVectorEraser(size_t f, size_t b) : front(f), back(b) {}
      template <class E>
      void visit(E& elem)
      {
        elem.erase( elem.begin()+front, elem.begin()+back );
      }
    };
    struct MultiVectorPrinter
    {
      std::ostream & s;
      size_t pos;
      std::string del;
      MultiVectorPrinter(std::ostream & _s, size_t p) : s(_s), pos(p), del("") {}
      template <class E>
      void visit(E& elem)
      {
        s << del << elem[pos];
        del = ", ";
      }
    };
  } // end empty namespace
#endif

  /**
     proxy object to give access

     @tparam T the tuple< vector<...> > type
   */
  template< typename T >
  struct MultiDataProxy;

  template< typename T >
  struct MultiDataProxy
  {
    typedef MultiDataProxy<typename remove_const<T>::type> MutableProxy;
    typedef MultiDataProxy<const typename remove_const<T>::type> ConstProxy;
    T & _vectors;

    int pos;
    std::string name;
    MultiDataProxy(T & v, size_t pos, std::string _n) :
      _vectors(v), pos(pos), name(_n) {}
    MultiDataProxy(ConstProxy & other) :
      _vectors(other._vectors), pos(other.pos), name(other.name) {}
    MultiDataProxy(MutableProxy & other) :
      _vectors(other._vectors), pos(other.pos), name(other.name) {}
    // compare
    bool operator == (const MultiDataProxy & other) const { return get<0>() == other.get<0>(); }
    bool operator != (const MultiDataProxy & other) const { return get<0>() != other.get<0>(); }
    bool operator  < (const MultiDataProxy & other) const { return get<0>()  < other.get<0>(); }
    bool operator  > (const MultiDataProxy & other) const { return get<0>()  > other.get<0>(); }
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
    template <size_t N>
    typename TypeTraits<typename tuple_element<N,T>::type>::ReferredType::reference
    get() {
      return Dune::get<N>(_vectors)[pos];
    }
    template <size_t N>
    typename TypeTraits<typename tuple_element<N,T>::type>::ReferredType::const_reference
    get() const {
      return Dune::get<N>(_vectors)[pos];
    }
  private:
    template<typename P>
    void assign(const P & other) {
#ifdef DEBUG_MULTIVEC
      std::cerr << "Assign " << name << "," << pos
                << "\n   from " << other.name << "," << other.pos << std::endl;
#endif
#warning TODO: Assign-TMP
      get<0>() = other.get<0>();
      get<1>() = other.get<1>();
      get<2>() = other.get<2>();
      get<3>() = other.get<3>();
    }
  };

  template< typename T >
  std::ostream& operator<< (std::ostream & s, const MultiDataProxy<T> & i)
  {
    s << "(";
    MultiVectorPrinter printer(s, i.pos);
    ForEachValue<T> forEach(i._vectors);
    forEach.apply(printer);
    s << ")";
    return s;
  }

  template< typename T >
  class MultiVectorIterator :
    public Dune::BidirectionalIteratorFacade< MultiVectorIterator<T>,
        MultiDataProxy<T> >
  {
    // friend class MultiVectorIterator<typename remove_const<C>::type, typename remove_const<T>::type >;
    // friend class TestIterator<const typename remove_const<C>::type, const typename remove_const<T>::type >;
  public:
    mutable MultiDataProxy<T> data;

    // constructors
    MultiVectorIterator(T & v, size_t n, std::string i) :
      data(v,n,i) {}
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
      assert(other.data._vectors == data._vectors);
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
      assert(other.data._vectors == data._vectors);
      return other.data.pos == data.pos;
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
    MultiDataProxy<T>& dereference() const
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

    typedef tuple<A&, B&, C&, D&> T;
    T _vectors;
    std::string _name;

  public:
    //! container interface typedefs
    //! \{

    /** \brief Type of the values stored by the container */
    typedef MultiDataProxy<T> value_type;

    /** \brief Type of the const values stored by the container */
    typedef MultiDataProxy<T> const_value_type;

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
    typedef MultiVectorIterator<T> iterator;
    typedef MultiVectorIterator<const T> const_iterator;
    //! \}

    MultiVector(A & a, B & b, C & c, D & d, std::string n = "?") :
      _vectors(tie(a,b,c,d))
    {
      assertSize();
      _name = n;
    }

    iterator begin()
    {
      return iterator(_vectors, 0, _name);
    }

    iterator end()
    {
      assertSize();
      return iterator(_vectors, get<0>().size(), _name);
    }

    void erase(iterator front, iterator back)
    {
      assertSize();

      MultiVectorEraser erase(front.pos(), back.pos());
      ForEachValue<T> forEach(_vectors);
      forEach.apply(erase);
    }
    size_t size() const {
      assertSize();
      return get<0>().size();
    }

    template <size_t N>
    typename tuple_element<N,T>::type &
    get() {
      return Dune::get<N>(_vectors);
    }
    template <size_t N>
    typename tuple_element<N,T>::type
    get() const {
      return Dune::get<N>(_vectors);
    }
  private:
    void assertSize() const {
      MultiVectorAssertSize check(get<0>().size());
      ForEachValue<const T> forEach(_vectors);
      forEach.apply(check);
    }
  };

} // end namespace Dune

#endif // DUNE_MULTIVECTOR_HH
