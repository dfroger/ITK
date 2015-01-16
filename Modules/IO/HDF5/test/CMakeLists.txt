# todo:
# itkHDF5ImageIOStreamingReadWriteTest.cxx
# itkHDF5ImageIOPastingFromMemoryTest.cxx
# itkHDF5ImageIOPastingFromFileTest.cxx
# itkHDF5ImageIOPastingInPlaceTest.cxx
# Plus: with streaming

itk_module_test()
set(ITKIOHDF5Tests
  itkHDF5ImageIOTest.cxx
  itkHDF5ImageIOStreamingReadWriteTest.cxx
  itkHDF5ImageIOPastingFromMemoryTest.cxx
  itkHDF5ImageIOPastingFromFileTest.cxx
)

CreateTestDriver(ITKIOHDF5  "${ITKIOHDF5-Test_LIBRARIES}" "${ITKIOHDF5Tests}")

itk_add_test(NAME itkHDF5ImageIOTest
  COMMAND ITKIOHDF5TestDriver itkHDF5ImageIOTest ${ITK_TEST_OUTPUT_DIR} )

itk_add_test(NAME itkHDF5ImageIOStreamingReadWriteTest
  COMMAND ITKIOHDF5TestDriver itkHDF5ImageIOStreamingReadWriteTest ${ITK_TEST_OUTPUT_DIR} )

itk_add_test(NAME itkHDF5ImageIOPastingFromMemoryTest
  COMMAND ITKIOHDF5TestDriver itkHDF5ImageIOPastingFromMemoryTest ${ITK_TEST_OUTPUT_DIR} )

itk_add_test(NAME itkHDF5ImageIOPastingFromFileTest
  COMMAND ITKIOHDF5TestDriver itkHDF5ImageIOPastingFromFileTest ${ITK_TEST_OUTPUT_DIR} )
