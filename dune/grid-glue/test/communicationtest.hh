// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef GRIDGLUE_COMMTEST_HH
#define GRIDGLUE_COMMTEST_HH

template <typename ctype, int dimw, bool forward>
class CheckGlobalCoordDataHandle :
  public Dune::GridGlue::CommDataHandle< CheckGlobalCoordDataHandle<ctype, dimw, forward>, Dune::FieldVector<ctype,dimw> >
{
  template<class RISType>
  bool hasSender (RISType& i) const
  {
    if (forward)
      return i.self();
    else
      return i.neighbor();
  }
  template<class RISType>
  bool hasRecipient (RISType& i) const
  {
    if (forward)
      return i.neighbor();
    else
      return i.self();
  }
public:
  template<class RISType>
  size_t size (RISType& i) const
  {
    assert(hasSender(i));
    if (forward)
      return i.geometry().corners();
    else
      return i.geometryOutside().corners();
  }

  template<class MessageBuffer, class EntityType, class RISType>
  void gather (MessageBuffer& buff, const EntityType& e, const RISType & i) const
  {
    assert(hasSender(i));
    if (forward)
      for (size_t n=0; n<size(i); n++)
        buff.write(i.geometry().corner(n));
    else
      for (size_t n=0; n<size(i); n++)
        buff.write(i.geometryOutside().corner(n));
  }

  template<class MessageBuffer, class EntityType, class RISType>
  void scatter (MessageBuffer& buff, const EntityType& e, const RISType & i, size_t n)
  {
    assert(hasRecipient(i));
    assert(n == size(i));
    for (size_t n=0; n<size(i); n++)
    {
      Dune::FieldVector<ctype,dimw> x;
      buff.read(x);
      if (forward)
        assert( (x - i.geometry().corner(n)).two_norm() < 1e-6 );
      else
        assert( (x - i.geometryOutside().corner(n)).two_norm() < 1e-6 );
    }
  }
};

template <class GlueType>
void testCommunication (const GlueType& glue)
{
  typedef typename GlueType::ctype ctype;
  enum { dimw = GlueType::dimworld };
  CheckGlobalCoordDataHandle<ctype, dimw, true> dh_forward;
  CheckGlobalCoordDataHandle<ctype, dimw, false> dh_backward;
  glue.communicate(dh_forward, Dune::All_All_Interface, Dune::ForwardCommunication);
  glue.communicate(dh_backward, Dune::All_All_Interface, Dune::BackwardCommunication);
}

#endif // GRIDGLUE_COMMTEST_HH
