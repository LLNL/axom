#ifndef SINA_RELATIONSHIP_HPP
#define SINA_RELATIONSHIP_HPP

/**
 * @file
 *
 * Contains the definition of the Relationship class.
 */

#include <string>

#include "conduit.hpp"

#include "sina/ID.hpp"

namespace sina {

/**
 * A Relationship consists of three parts: a subject, an object, and a predicate. It describes
 * a relationship between two Records (and/or Record inheritors, e.g. Run).
 * The subject and object must be IDs referring to valid records, while the predicate may be
 * any string.
 *
 * In describing the connection between objects, a Relationship is read as
 * "<subject> <predicate> <object>". For example, in the relationship
 * "Alice knows Bob", "Alice" is the subject, "knows" is the predicate, and
 * "Bob" is the object. For further examples:
 *
 * - Task_22 contains Run_1024
 * - msub_1_1 describes out_j_1_1
 * - Carlos sends an email to Dani
 * - local_task_12 runs before local_run_14
 *
 * Note that Relationships are described in the active voice. **Avoiding the passive voice
 * in predicates is recommended**, as this keeps the "direction" of the relationship constant.
 * An example of a passively-voiced Relationship is "Dani is emailed by Carlos". Instead,
 * this should be phrased as "Carlos emails Dani".
 *
 * If assembling Relationships programatically, it may be useful to reference the
 * ID documentation.
 *
 * \code
 * sina::ID task22{"Task_22", sina::IDType::Global};
 * sina::ID run1024{"Run_1024", sina::IDType::Global};
 * sina::Relationship myRelationship{task22, "contains", run1024};
 * std::cout << myRelationship.toNode().to_json() << std::endl;
 * \endcode
 *
 * This would output:
 * \code{.json}
 * {"object":"Run_1024","predicate":"contains","subject":"Task_22"}
 * \endcode
 *
 * As with any other Sina ID, the subject or object may be either local (uniquely refer to one object
 * in a Sina file) or global (uniquely refer to one object in a database). Local IDs are replaced with
 * global ones upon ingestion; all Relationships referring to that Local ID (as well as the Record possessing
 * that ID) will be updated to use the same global ID.
 *
 * \code
 * sina::ID myLocalID{"my_local_run", sina::IDType::Local};
 * std::unique_ptr<sina::Record> myRun{new sina::Run{myLocalID, "My Sim Code", "1.2.3", "jdoe"}};
 * sina::Relationship myRelationship{task22, "contains", myLocalID};
 * \endcode
 *
 * In the above code, "my_local_run" would be replaced by a global ID on ingestion. If this new global ID was,
 * for example, "5Aed-BCds-23G1", then "my_local_run" would automatically be replaced by "5Aed-BCds-23G1" in both
 * the Record and Relationship entries.
 */
class Relationship {
public:
    /**
     * Create a new relationship.
     *
     * @param subject the subject of the relationship
     * @param predicate the predicate describing the relationship from the
     * subject to the object
     * @param object the object of the relationship
     */
    Relationship(ID subject, std::string predicate, ID object);

    /**
     * Create a Relationship object from its representation as a conduit Node.
     *
     * @param asNode the relationship as a Node
     */
    explicit Relationship(conduit::Node const &asNode);

    /**
     * Get the subject.
     *
     * @return the subject
     */
    ID const &getSubject() const noexcept {
        return subject.getID();
    }

    /**
     * Get the object.
     *
     * @return the object
     */
    ID const &getObject() const noexcept {
        return object.getID();
    }

    /**
     * Get the predicate.
     *
     * @return the predicate
     */
    std::string const &getPredicate() const noexcept {
        return predicate;
    }

    /**
     * Convert this Relationship to its Node representation.
     *
     * @return this relationship as a conduit Node
     */
    conduit::Node toNode() const;

private:
    internal::IDField subject;
    internal::IDField object;
    std::string predicate;
};

}

#endif //SINA_RELATIONSHIP_HPP
