/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef STLREADER_HPP_
#define STLREADER_HPP_

#include <string>
#include <vector>

#include "axom/Macros.hpp"
#include "mint/UnstructuredMesh.hpp"

namespace axom {  
namespace quest 
{

class STLReader
{
    static const std::size_t  BINARY_HEADER_SIZE = 80; // bytes
    static const std::size_t  BINARY_TRI_SIZE = 50;    // bytes

public:

    /*!
     ***************************************************************************
     * \brief Constructor.
     ***************************************************************************
     */
    STLReader();

    /*!
     ***************************************************************************
     * \brief Destructor.
     ***************************************************************************
     */
    virtual ~STLReader();

    /*!
     ***************************************************************************
     * \brief Sets the name of the file to read.
     * \param [in] fileName the name of the file to read.
     ***************************************************************************
     */
    void setFileName( const std::string& fileName ) { m_fileName = fileName; };

    /*!
     ***************************************************************************
     * \brief Clears all internal data-structures
     ***************************************************************************
     */
    void clear();

    /*!
     ***************************************************************************
     * \brief Reads in the surface mesh from an STL file.
     * \pre m_fileName != ""
     ***************************************************************************
     */
    virtual void read();

    /*!
     ***************************************************************************
     * \brief Stores the STL data in the supplied unstructured mesh object.
     * \param [in,out] mesh pointer to the unstructured mesh.
     * \pre mesh != AXOM_NULLPTR.
     ***************************************************************************
     */
    void getMesh( axom::mint::UnstructuredMesh< MINT_TRIANGLE >* mesh );


private:
    /**
     * \brief A predicate to check if the file is in ascii format
     */
    bool isAsciiFormat() const;

    /**
     * \brief Reads an ascii encoded STL file into memory
     * \note The filename should be set with STLReader::setFileName()
     */
    void readAsciiSTL();

    /**
     * \brief Reads a binary encoded STL file into memory
     * \note The filename should be set with STLReader::setFileName()
     */
    void readBinarySTL();

protected:
    std::string m_fileName;

    int m_num_nodes;
    int m_num_faces;

    std::vector<double> m_nodes;

private:

    DISABLE_COPY_AND_ASSIGNMENT(STLReader);
    DISABLE_MOVE_AND_ASSIGNMENT(STLReader);
};

} // end namespace quest 
} // end namespace axom 

#endif /* STLREADER_HPP_ */
