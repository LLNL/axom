#ifndef SINA_FILE_HPP
#define SINA_FILE_HPP

/// @file

#include <string>
#include <vector>

#include "conduit.hpp"

namespace sina {
/**
 * A File tracks the location (URI) and mimetype of a file on the file system, plus any tags.
 * In the Sina schema, a File always belongs to a Record or one of Record's inheriting types.
 *
 * Every File must have a URI, while mimetype and tags are optional.
 *
 * \code
 * sina::File myFile{"/path/to/file.png"};
 * myFile.setMimeType("image/png");
 * sina::File myOtherFile{"/path/to/other/file.txt"};
 * myOtherFile.setTags({"these","are","tags"});
 * myRecord->add(myFile);
 * myRecord->add(myOtherFile);
 * \endcode
 */
class File {
public:
    /**
     * Construct a new File.
     *
     * @param uri the location of the file
     */
    explicit File(std::string uri);

    /**
     * Construct a new File.
     *
     * @param uri the uri for a file
     * @param asNode the Node representation of the file's additional info
     */
    File(std::string uri, conduit::Node const &asNode);

    /**
     * Get the File's URI.
     *
     * @return the URI
     */
    std::string const &getUri() const noexcept {
        return uri;
    }

    /**
     * Get the File's MIME type.
     *
     * @return the MIME type
     */
    std::string const &getMimeType() const noexcept {
        return mimeType;
    }

    /**
     * Get the File's tags.
     *
     * @return the tags
     */
    std::vector<std::string> const &getTags() const noexcept {
        return tags;
    }

    /**
     * Set the File's MIME type.
     *
     * @param mimeType the MIME type
     */
    void setMimeType(std::string mimeType);

    /**
     * Set the File's tags.
     *
     * @param tags the File's tags
     */
    void setTags(std::vector<std::string> tags);

    /**
     * Convert this File to its conduit Node representation.
     *
     * @return the File in its Node representation
     */
    conduit::Node toNode() const;

private:
    std::string uri;
    std::string mimeType;
    std::vector<std::string> tags;
};

}

#endif //SINA_FILE_HPP
