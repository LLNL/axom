/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file STLReader.hpp
 *******************************************************************************
 */

#ifndef STLREADER_HPP_
#define STLREADER_HPP_

#include <string>
#include <vector>

#include "common/ATKMacros.hpp"
#include "quest/UnstructuredMesh.hpp"

namespace quest
{

class STLReader
{
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
     * \pre mesh != ATK_NULLPTR.
     ***************************************************************************
     */
    void getMesh( meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE >* mesh );


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
};

} /* namespace quest */

#endif /* STLREADER_HPP_ */
