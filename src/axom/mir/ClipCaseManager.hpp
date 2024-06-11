#define ST_TRI  0
#define ST_QUAD 1
#define ST_TET  2
#define ST_PYR  3
#define ST_WDG  4
#define ST_HEX  5
#define ST_MAX  6

template <typename ContainerType, MemorySpace SPACE>
struct ClipCaseBase
{
  using IntContainerType = ContainerType<int, 1, SPACE>;
  using UInt8ContainerType = ContainerType<uint8, 1, SPACE>;

  // Q: Do I need to explicitly provide a constructor to get it marked as AXOM_HOST_DEVICE?

  size_t size() const { return m_shapes.size(); }

  size_t shapesForCase(size_t index) const
  {
    return m_shapes[index];
  }

  UInt8ContainerType getShape(size_t caseId, size_t shapeId) const
  {
    assert(caseId < m_shapes.size());
    assert(shapeId < shapesForCase(caseId));

    const uint8 *shapeStart = m_table.data() + m_offsets[caseId];
    size_t shapeLen = 0;
    for(int i = 0; i < shapeId; i++)
    {
      shapeLen = advance(*shapeStart);
      shapeStart += shapeLen;
    }
    shapeLen = advance(*shapeStart);
    return UInt8ContainerType(shapeStart, shapeLen);
  }

  size_t advance(uint8 shape) const
  {
    if(shape == ST_TRI)
      retval = 2 + 3;
    else if(shape == ST_QUAD)
      retval = 2 + 4;
    else if(shape == ST_TET)
      retval = 2 + 4;
    else if(shape == ST_PYR)
      retval = 2 + 5;
    else if(shape == ST_WDG)
      retval = 2 + 6;
    else if(shape == ST_HEX)
      retval = 2 + 8;
    return retval;
  }
  
  IntContainerType   m_shapes;
  IntContainerType   m_offsets;
  UInt8ContainerType m_table;
};

template <MemorySpace SPACE>
struct ClipCase : public ClipCaseBase<axom::Array, SPACE>
{
  using ClipCaseView = ClipCaseBase<axom::ArrayView, SPACE>

  void load(size_t n, const int *shapes, const int *offsets, const uint8 *table, size_t tableLen)
  {
    m_shapes = IntContainerType(shapes, n);
    m_offsets = IntContainerType(offsets, n);
    m_table = UInt8ContainerType(table, tableLen);
  }

  ClipCaseView view() const
  {
    ClipCaseView v;
    v.m_shapes = m_shapes.view();
    v.m_offsets = m_offsets.view();
    v.m_table = m_table.view(); 
    return v;
  }

};


class ClipCaseManager
{
public:
  ClipCaseManager()
  {
    for(int i = 0; i < ST_MAX; i++)
      m_clipCases[i] = ClipCase<SPACE>();
  }

  const ClipCase &operator[](size_t shapeId)
  {
    assert(shapeID < ST_MAX);
    if(m_clipCases[shapeId].size() == 0)
    {
      load(shapeId);
    }
    return m_clipCases[shapeId];
  }

private:
  void load(size_t shapeId)
  {
    if(shapeId == ST_TRI)
      m_clipCases[shapeId].load(numClipCasesTri, numClipShapesTri, startClipShapesTri, clipShapesTri);
    else if(shapeId == ST_QUAD)
      m_clipCases[shapeId].load(numClipCasesQuad, numClipShapesQuad, startClipShapesQuad, clipShapesQuad);
    else if(shapeId == ST_TET)
      m_clipCases[shapeId].load(numClipCasesTet, numClipShapesTet, startClipShapesTet, clipShapesTet);
    else if(shapeId == ST_PYR)
      m_clipCases[shapeId].load(numClipCasesPyr, numClipShapesPyr, startClipShapesPyr, clipShapesPyr);
    else if(shapeId == ST_WDG)
      m_clipCases[shapeId].load(numClipCasesWdg, numClipShapesWdg, startClipShapesWdg, clipShapesWdg);
    else if(shapeId == ST_HEX)
      m_clipCases[shapeId].load(numClipCasesHex, numClipShapesHex, startClipShapesHex, clipShapesHex);
  }

  ClipCase<SPACE> m_clipCases[ST_MAX];
};
