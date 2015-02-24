% *************************************************************************
% function demo_hashtable
%
% Purpose: shows the usage of almost all mHashTable methods. All mHashTable
%          methods follow closely Java, see
%          http://java.sun.com/j2se/1.4.2/docs/api/java/util/Hashtable.html
%
% *************************************************************************
function demo_hashtable

clc
N       = 2000;

% create mHashTable and Java Hashtable objects
my_h    = mhashtable;
ja_h    = java.util.Hashtable;

disp(' ');
disp('-put dummy data into mHashTable and get elapsed time (not reliable)')
tic
for i = 1:N, my_h.put(i*pi, i); end
toc


disp(' ');
disp('-put dummy data into Java Hashtable and get elapsed time (not reliable)')
tic
for i = 1:N, ja_h.put(i*pi, i); end
toc


disp(' ');
disp('-clears this hashtable so that it contains no keys')
my_h.clear


disp(' ');
disp('-put different kinds of data to the REFERENCE and dispaly contents')
s = local_PutDummyData(my_h);
my_h.toString

disp(' ');
disp('-get size')
my_h.size

disp(' ');
disp('-try to use a struct as the key - not allowed!!! success == 0');
success = my_h.put(s, s)


disp(' ');
disp('-check whether hashtable contains a value');
my_h.contains('short string')


disp(' ');
disp('-check whether hashtable contains a key');
my_h.containsKey(uint8(45))


disp(' ');
disp('-return all values');
my_h.elements


disp(' ');
disp('-check whether two mHashTables equals');
my_h.equals(my_h)


disp(' ');
disp('-get non-stored the value from hashtable (returns empty)');
non_existing_value = my_h.get('non-existing double array')


disp(' ');
disp('-get a stored value');
double_array = my_h.get('double array')


disp(' ');
disp('-get hash index of this object (actually of no use)');
my_h.hashcode


disp(' ');
disp('-check whether hashtable is empty');
my_h.isempty


disp(' ');
disp('-get all keys');
my_h.keys


disp(' ');
disp('-remove the element corresponding to the keys ''cloned copy'' and 45');
my_h.remove('cloned copy');
my_h.remove(45);
my_h


disp(' ');
disp('-use a non-java method ''getElement'' and return the mHashTable intern representation');
my_h.getElement

disp(' ');
disp('-finally delete "my_h" and free the memory');
my_h.delete     % Syntax delete(my_h) is also possible

disp(' ');
disp('-check whether "my_h" is still a valid mHashTable obect');
my_h.isMHashTable

disp(' ');
disp('-"my_h" does not exist in memory any more, see what happens if it is nevertheless used');
some_value = my_h.get(1)


% *************************************************************************
% function s = local_PutDummyData(my_h)
%
% Purpose: put different types of data into the hashtable, return a dummy
%          struct
%
% *************************************************************************
function s = local_PutDummyData(my_h)

% prepare dummy data
N       = 10;
cl      = cell(N, 1);
s(1).w  = {'string 1', 'string 2'};
s(1).d  = [1,2,3,4,5,6];
s(2).w  = {434, 744};
s(2).d  = 'this is a string';
for i = 1:N, cl{i} = sprintf('%f', i); end

my_h.put(cl{1}, 'short string');
my_h.put(cl{2}, s(1));
my_h.put(cl{3}, {'cell array string', 23, s});
my_h.put(cl{4}, s);
my_h.put('"hallo world" key with "pi" value', pi);
my_h.put(1, 'another short string');
my_h.put(exp(1), 'very long string - so long that only [1x67 char] would be displayed');
my_h.put('cloned copy', my_h.clone);
my_h.put('double array', 1e-9*rand(1,4));
my_h.put(uint8(45), 'the key is an uint8(45)');
my_h.put('small array of uint16', uint16([1:2]));
my_h.put('large array of uint32', uint32([1:20]));