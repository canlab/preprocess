mHashTable Release 1.1, Dec 2005
by Mikhail Poda-Razumtsev (Munich, Germany)

mHashTable is free software under the GNU General Public License,
http://www.gnu.org/copyleft/gpl.html

Requirements: MATLAB Release 6.5 - 7.1


WHAT IS A HASH TABLE?

Hash table, or a hash map, is a universal data structure that associates keys with values.
The primary operation it supports efficiently is a lookup: given a key
(e.g. a person's name), find the corresponding value (e.g. that person's telephone number).
It works by transforming the key using a hash function into a hash,
a number that the hash table uses to locate the desired value. Usually both keys and values
are allowed to be of arbitrary data type.


MATLAB HASH TABLE

When building a medium or large-size Matlab-based  application there is often a need for the
universal data container. The data can be of any type - struct, cell or double.
This data should be accessed through keys which can by also of any type.

The easiest way to implement such container with Matlab would be using a struct array
or a cell array because both can contain any matlab data type.
But both approaches have serious limitations:

1. Struct array with the keys as its fields is suitable only if the keys are Matlab variables
2. Cell array with the keys as indices is suitable only if the keys is a series of integers
3. None approach is suitable if the keys are of both types - strings and numbers
4. None approach creates a REFERENCE, because matlab does not support references (unlike C)

Other programming languages support such universal data containers.
For Java it is the Java Collections Framework. See a very good introduction:
http://www.digilife.be/quickreferences/PT/Java%20Collections%20Framework.pdf
Besides lists and sets it contains a hashtable class.

Java Hashtable can be used in Matlab. Its sintax is very easy and it can be used by inexperienced 
programmers. And it has a big advantage - it is a REFERENCE (all Java objects are references)!
But the problem using java Hashtable in matlab is that it supports not all matlab data types
but only strings and numbers. Matlab structures and cell arrays can not be stored in the Java Hashtable.


@mHashTable  IMPLEMENTATION

See the example "demo_hashtable.m"

The mHashTable class is designed to be for the matlab what Hashtable is for the Java.
It is a reference, it has Java syntax and it accepts values of any type.
All java Hashtable methods are implemented for the mHashTable. For the Java Hashtable
methods description see: http://java.sun.com/j2se/1.4.2/docs/api/java/util/Hashtable.html

The reference is achieved by storing all data in a persistent variable in the mHashTable
private function called HashtableStore, the mHashTable objects contain only the index 
to an element in this persistent variable.

The true Hashtable has constant access time, independent of the Hashtable size. 
This feature is difficult/impossible to implement/would be very slow in matlab.
mHashTable uses linear search to retrive the keys, thus the access time 
grows proporitonally to the Hashtable size.

The key are allowed to be numeric or string, but this is not a real limitation because other
types are usually of no interest.

If trying to get the data from non-present key then an empty array is returned.

mHashTable class is implemented entirely in Matlab with the aim of being so fast as matlab would allow,
but it is still 2-5 times slower then corresponding Java Hashtable for moderate sizes.


IMPORTANT

There is no automatic garbage collection. Please use the nice "delete" method to free the 
hashtable objects which you are not using any more. Usually it is not a problem if 
you are not doing it, but in some cases it would speed up the access time 
to other hashtables (if you have some).