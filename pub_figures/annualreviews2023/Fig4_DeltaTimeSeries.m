   l   i   s   t       =       [   "   D   v   i   n   a   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   a   m   p   u   t   r   a   "   ,   "   M   a   c   K   e   n   z   i   e   "   ,   "   M   e   k   o   n   g   "   ,   "   M   i   s   s   i   s   s   i   p   p   i   "   ,   "   N   i   g   e   r   "   ,   "   N   i   l   e   "   ,   "   P   a   r   a   n   a   "   ,   "   W   a   i   p   a   o   a   "   ,   "   Y   e   l   l   o   w   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   m   a   p   u   t   r   a   "   ,   "   S   o   n   g       H   o   n   g   "   ,   "   P   o   "   ,   "   R   h   i   n   e   -   M   e   u   s   e   "   ,   "   T   a   m   a   "   ,   "   T   r   i   n   i   t   y   "   ]   ;   
   
   
   t   i   l   e   d   l   a   y   o   u   t   (   '   f   l   o   w   '   )   
   n   e   x   t   t   i   l   e   
   
   %   p   a   s   t       S   L   
   d       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   a   n   u   _   r   s   l   \   d   e   l   t   a   _   r   s   l   _   a   n   u   h   r   _   7   1   p   2   3   0   .   d   a   t   '   )   ;   
   s   l   _   t   i   m   e       =       d   {   :   ,   1   }   ;   
   s   l   _   d   a   t   a       =       d   {   :   ,   [   2   :   2   :   e   n   d   ]   }   ;   
   
   %   f   u   t   u   r   e       S   L   
   l   o   a   d   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   d   e   l   t   a   -   t   s   .   m   a   t   '   )   
   
   t       =       [   s   l   _   t   i   m   e   (   2   :   e   n   d   )       ;       d   o   u   b   l   e   (   t   i   m   e   (   2   :   e   n   d   )   )   .   /   1   0   0   0   -   2   ]   ;   
   s   l   r   a   t   e       =       [   d   i   f   f   (   s   l   _   d   a   t   a   )   .   /   d   i   f   f   (   s   l   _   t   i   m   e   )   ;       d   i   f   f   (   S   S   P   2   4   5   _   5   0   (   2   :   e   n   d   ,   4   :   e   n   d   )   '   )   .   /   1   0   ]   ;   
   
   
   
   p   l   o   t   (   t   ,   s   l   r   a   t   e   )   
   y   l   a   b   e   l   (   '   R   S   L   R       (   m   m   /   y   r   )   '   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   L   i   m   '   ,   [   -   2   0       2   0   ]   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   
   %   s   t   a   n   l   e   y       w   a   r   n   e   
   n   e   x   t   t   i   l   e   
   d       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   s   t   a   n   l   e   y   _   w   a   r   n   e   .   x   l   s   x   '   )   ;   
   h   i   s   t   o   g   r   a   m   (   -   d   .   A   g   e   O   f   _   n   e   a   r   _   B   a   s   a   l   H   o   l   o   c   e   n   e   _   y   e   a   r   s   _   .   /   1   0   0   0   ,   [   -   2   5   :   1   :   5   ]   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   y   l   a   b   e   l   (   '   N   u   m   b   e   r       o   f       d   e   l   t   a   s   '   )   
   b   o   x       o   n   
   
   n   e   x   t   t   i   l   e   
   %   q   u   o   c   k       d   a   t   a   
   t       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   q   u   o   c   k   .   x   l   s   x   '   )   ;   
   [   l   ,   i   d   x   ]       =       u   n   i   q   u   e   (   t   .   B   a   s   i   n   I   D   2   )   ;   
   f   o   r       i   i   =   1   :   l   e   n   g   t   h   (   l   )   ,   
                   h   o   l   d       o   n   
                   d       =       t   .   D   i   s   t   a   n   c   e   _   k   m   (   t   .   B   a   s   i   n   I   D   2   =   =   l   (   i   i   )   )   ;   
                   p   l   o   t   (   -   t   .   T   i   m   e   _   k   y   a   B   P   (   t   .   B   a   s   i   n   I   D   2   =   =   l   (   i   i   )   )   ,   d   .   /   d   (   e   n   d   )   .   *   1   0   0   )   
   e   n   d                   
   l   e   g   e   n   d   (   t   .   D   e   l   t   a   (   i   d   x   )   ,   '   L   o   c   a   t   i   o   n   '   ,   '   S   o   u   t   h   W   e   s   t   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   2   5       5   ]   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   y   l   a   b   e   l   (   '   D   e   l   t   a       L   e   n   g   t   h       c   o   m   p   a   r   e   d       t   o       m   o   d   e   r   n       (   %   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   b   o   x       o   n   
   
   
   s   a   v   e   a   s   (   g   c   f   ,   '   F   i   g   4   _   D   e   l   t   a   T   i   m   e   S   e   r   i   e   s   2   .   p   n   g   '   )   
   
   %   %       a   l   t   e   r   n   a   t   i   v   e   
   c   l   r   
   
   l   i   s   t       =       [   "   D   v   i   n   a   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   a   m   p   u   t   r   a   "   ,   "   M   a   c   K   e   n   z   i   e   "   ,   "   M   e   k   o   n   g   "   ,   "   M   i   s   s   i   s   s   i   p   p   i   "   ,   "   N   i   g   e   r   "   ,   "   N   i   l   e   "   ,   "   P   a   r   a   n   a   "   ,   "   W   a   i   p   a   o   a   "   ,   "   Y   e   l   l   o   w   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   m   a   p   u   t   r   a   "   ,   "   S   o   n   g       H   o   n   g   "   ,   "   P   o   "   ,   "   R   h   i   n   e   -   M   e   u   s   e   "   ,   "   T   a   m   a   "   ,   "   T   r   i   n   i   t   y   "   ]   ;   
   
   f       =       '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   '   ;   
   
   t   i   l   e   d   l   a   y   o   u   t   (   '   f   l   o   w   '   )   
   n   e   x   t   t   i   l   e   
   
   %   p   a   s   t       S   L   
   %   d       =       r   e   a   d   t   a   b   l   e   (   [   f       '   a   n   u   _   r   s   l   \   d   e   l   t   a   _   r   s   l   _   a   n   u   h   r   _   7   1   p   2   3   0   .   d   a   t   '   ]   )   ;   
   d       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   i   c   e   6   g   \   d   e   l   t   a   _   r   s   l   _   i   c   e   6   g   _   v   m   5   a   2   .   d   a   t   '   )   ;   
   s   l   _   t   i   m   e       =       d   {   :   ,   1   }   ;   
   s   l   _   d   a   t   a       =       d   {   :   ,   [   2   :   2   :   e   n   d   ]   }   ;   
   
   f   i   d       =       f   o   p   e   n   (   [   f       '   i   c   e   6   g   \   g   m   s   l   _   d   e   f   .   p   l   o   t   '   ]   )   ;       
   g   m   s   l   _   d   e   f       =       t   e   x   t   s   c   a   n   (   f   i   d   ,   '   %   f       %   f   '   )   ;       
   g   m   s   l   _   d   e   f       =       [   g   m   s   l   _   d   e   f   {   :   }   ]   ;   
   f   c   l   o   s   e   (   f   i   d   )   ;   
   
   f   i   d       =       f   o   p   e   n   (   [   f       '   i   c   e   6   g   \   g   m   s   l   _   m   w   .   p   l   o   t   '   ]   )   ;       
   g   m   s   l   _   m   w       =       t   e   x   t   s   c   a   n   (   f   i   d   ,   '   %   f       %   f   '   )   ;       
   g   m   s   l   _   m   w       =       [   g   m   s   l   _   m   w   {   :   }   ]   ;   
   f   c   l   o   s   e   (   f   i   d   )   ;   
   
   f   i   d       =       f   o   p   e   n   (   [   f       '   i   c   e   6   g   \   g   m   s   l   _   n   e   t   .   p   l   o   t   '   ]   )   ;       
   g   m   s   l   _   n   e   t       =       t   e   x   t   s   c   a   n   (   f   i   d   ,   '   %   f       %   f   '   )   ;       
   g   m   s   l   _   n   e   t       =       [   g   m   s   l   _   n   e   t   {   :   }   ]   ;   
   f   c   l   o   s   e   (   f   i   d   )   ;   
   
   %   f   u   t   u   r   e       S   L   
   l   o   a   d   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   d   e   l   t   a   -   t   s   .   m   a   t   '   )   
   
   t       =       s   l   _   t   i   m   e   (   2   :   e   n   d   )   ;   
   
   s   l   r   a   t   e       =       d   i   f   f   (   s   l   _   d   a   t   a   )   .   /   d   i   f   f   (   s   l   _   t   i   m   e   )   ;   
   p   l   o   t   (   t   ,   s   l   r   a   t   e   (   :   ,   [   1   1   :   1   6   ]   )   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   ;       %   ,   '   C   o   l   o   r   '   ,   [   0   .   5       0   .   5       0   .   5   ]   )   
   
   h   o   l   d       o   n   
   
   g   m   s   l       =       d   i   f   f   (   g   m   s   l   _   n   e   t   (   :   ,   2   )   )   .   /   d   i   f   f   (   g   m   s   l   _   n   e   t   (   :   ,   1   )   )   ;   
   p   l   o   t   (   g   m   s   l   _   n   e   t   (   2   :   e   n   d   ,   1   )   ,   g   m   s   l   ,   '   k   '   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   
   m   i       =       [   -   7   ,   -   8   ,   -   7   .   5   ,   -   5   ,   -   7   .   5   ,   -   2   .   5   ]   ;       %   t   i   m   e       i   n       k   a       B   P       o   f       m   i   n   i   m   u   m       i   n       e   x   t   e   n   t       f   r   o   m       q   u   o   c   k   
   [   ~   ,   i   s   a   ]       =       i   s   m   e   m   b   e   r   (   m   i   ,   t   )   ;   
   a       =       f   l   i   p   u   d   (   g   e   t   (   g   c   a   ,   '   C   h   i   l   d   r   e   n   '   )   )   ;       a       =       g   e   t   (   a   ,   '   C   o   l   o   r   '   )   ;   
   f   o   r       i   i   =   1   :   6   ,   
   s   c   a   t   t   e   r   (   t   (   i   s   a   (   i   i   )   )   ,   s   l   r   a   t   e   (   i   s   a   (   i   i   )   ,   1   0   +   i   i   )   ,   5   0   ,   '   M   a   r   k   e   r   F   a   c   e   C   o   l   o   r   '   ,   c   e   l   l   2   m   a   t   (   a   (   i   i   )   )   ,   '   M   a   r   k   e   r   E   d   g   e   C   o   l   o   r   '   ,   '   k   '   )   
   e   n   d   
   
   l   e   g   e   n   d   (   l   i   s   t   (   1   1   :   1   6   )   ,   '   L   o   c   a   t   i   o   n   '   ,   '   S   o   u   t   h   W   e   s   t   '   )   
   
   
   
   s   l   _   f   u   t   u   r   e       =       d   i   f   f   (   [   S   S   P   1   2   6   _   5   0   (   :   ,   4   :   e   n   d   )   ;       S   S   P   2   4   5   _   5   0   (   :   ,   4   :   e   n   d   )   ;   S   S   P   5   8   5   _   5   0   (   :   ,   4   :   e   n   d   )   ]   '   )   .   /   1   0   ;   
   t       =       [   0   ;       d   o   u   b   l   e   (   t   i   m   e   (   1   :   e   n   d   -   1   )   )   .   /   1   0   0   0   -   2   ]   ;   
   %   p   l   o   t   (   t   ,   [   r   e   p   m   a   t   (   s   l   r   a   t   e   (   e   n   d   ,   :   )   ,   1   ,   3   )   ;   s   l   _   f   u   t   u   r   e   (   :   ,   [   1   :   1   0   ,   1   2   :   2   1   ,   2   3   :   3   2   ]   )   ]   ,   '   C   o   l   o   r   '   ,   [   0   .   5       0   .   5       0   .   5   ]   )   
   p   l   o   t   (   t   ,   [   r   e   p   m   a   t   (   g   m   s   l   (   e   n   d   )   ,   1   ,   3   )   ;   s   l   _   f   u   t   u   r   e   (   :   ,   [   1   1       2   2       3   3   ]   )   ]   ,   '   k   '   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   
   y   l   a   b   e   l   (   '   R   S   L   R       (   m   m   /   y   r   )   '   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   1       0   .   5   ]   )   ;   
   
   n   e   x   t   t   i   l   e   
   %   q   u   o   c   k       d   a   t   a   
   t       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   q   u   o   c   k   .   x   l   s   x   '   )   ;   
   [   l   ,   i   d   x   ]       =       u   n   i   q   u   e   (   t   .   B   a   s   i   n   I   D   2   )   ;   
   f   o   r       i   i   =   1   :   l   e   n   g   t   h   (   l   )   ,   
                   h   o   l   d       o   n   
                   d       =       t   .   n   o   n   d   i   m   e   s   i   o   n   a   l   i   z   e   d   D   i   s   t   a   n   c   e   _   d   i   s   t   a   n   c   e   A   t   T   _   d   i   s   t   a   n   c   e   T   0   _   (   t   .   B   a   s   i   n   I   D   2   =   =   l   (   i   i   )   )   ;   
                   t   t       =       -   t   .   T   i   m   e   _   k   y   r   B   P   (   t   .   B   a   s   i   n   I   D   2   =   =   l   (   i   i   )   )   ;   
                   p   l   o   t   (   t   t   ,   d   *   1   0   0   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   ,   '   C   o   l   o   r   '   ,   a   {   i   i   }   )   ;   
                   %   [   ~   ,   j   j   ]       =       m   i   n   (   d   )   ;   
                   %   s   c   a   t   t   e   r   (   t   t   (   j   j   )   ,   d   (   j   j   )   .   *   1   0   0   ,   '   M   a   r   k   e   r   F   a   c   e   C   o   l   o   r   '   ,   g   e   t   (   a   ,   '   C   o   l   o   r   '   )   )   
                   %   t   t   (   j   j   )   
                   
   e   n   d                   
   l   e   g   e   n   d   (   t   .   D   e   l   t   a   (   i   d   x   )   ,   '   L   o   c   a   t   i   o   n   '   ,   '   S   o   u   t   h   W   e   s   t   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   1       0   .   5   ]   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   y   l   a   b   e   l   (   '   D   e   l   t   a       L   e   n   g   t   h       c   o   m   p   a   r   e   d       t   o       m   o   d   e   r   n       (   %   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   
   
   %   s   t   a   n   l   e   y       w   a   r   n   e   
   n   e   x   t   t   i   l   e   
   d       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   s   t   a   n   l   e   y   _   w   a   r   n   e   .   x   l   s   x   '   )   ;   
   h   i   s   t   o   g   r   a   m   (   -   d   .   C   a   l   i   b   r   a   t   e   d   A   g   e   _   c   a   l   Y   r   B   P   _   1   9   5   0   _   (   [   1   :   3   0   ,   3   2   :   e   n   d   ]   )   .   /   1   0   0   0   ,   [   -   1   0   :   0   .   5   :   0   .   5   ]   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   c   a   l   i   b   r   a   t   e   d       k   a       B   P   )   '   )   
   y   l   a   b   e   l   (   '   N   u   m   b   e   r       o   f       d   e   l   t   a   s   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   1       0   .   5   ]   )   ;   
   
   
   s   a   v   e   a   s   (   g   c   f   ,   '   F   i   g   4   _   D   e   l   t   a   T   i   m   e   S   e   r   i   e   s   .   s   v   g   '   )   
   
   
   
   %   %       f   o   r       T   O   R   
   
   c   l   r   
   
   l   i   s   t       =       [   "   D   v   i   n   a   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   a   m   p   u   t   r   a   "   ,   "   M   a   c   K   e   n   z   i   e   "   ,   "   M   e   k   o   n   g   "   ,   "   M   i   s   s   i   s   s   i   p   p   i   "   ,   "   N   i   g   e   r   "   ,   "   N   i   l   e   "   ,   "   P   a   r   a   n   a   "   ,   "   W   a   i   p   a   o   a   "   ,   "   Y   e   l   l   o   w   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   m   a   p   u   t   r   a   "   ,   "   S   o   n   g       H   o   n   g   "   ,   "   P   o   "   ,   "   R   h   i   n   e   -   M   e   u   s   e   "   ,   "   T   a   m   a   "   ,   "   T   r   i   n   i   t   y   "   ]   ;   
   i   d   x       =       [   1   1   :   1   6   ]   ;   
   t   i   l   e   d   l   a   y   o   u   t   (   '   f   l   o   w   '   )   
   
   
   %   p   a   s   t       S   L   
   d       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   a   n   u   _   r   s   l   \   d   e   l   t   a   _   r   s   l   _   a   n   u   h   r   _   7   1   p   2   3   0   .   d   a   t   '   )   ;   
   s   l   _   t   i   m   e       =       d   {   :   ,   1   }   ;   
   s   l   _   d   a   t   a       =       d   {   :   ,   [   2   :   2   :   e   n   d   ]   }   ;   
   t       =       s   l   _   t   i   m   e   (   2   :   e   n   d   )   ;   
   
   d   2       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   i   c   e   6   g   \   d   e   l   t   a   _   r   s   l   _   i   c   e   6   g   _   v   m   5   a   2   .   d   a   t   '   )   ;   
   s   l   _   t   i   m   e   2       =       d   2   {   :   ,   1   }   ;   
   s   l   _   d   a   t   a   2       =       d   2   {   :   ,   [   2   :   2   :   e   n   d   ]   }   ;   
   t   2       =       s   l   _   t   i   m   e   2   (   2   :   e   n   d   )   ;   
   
   n   e   x   t   t   i   l   e   
   
   p   l   o   t   (   s   l   _   t   i   m   e   ,   s   l   _   d   a   t   a   (   :   ,   i   d   x   )   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   h   o   l   d       o   n   
   t   i   t   l   e   (   '   A   N   U   '   )   
   y   l   a   b   e   l   (   '   R   S   L       (   m   )   '   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   2       0   ]   ,   '   Y   L   i   m   '   ,   [   -   6   0       1   0   ]   )   
   
   
   
   n   e   x   t   t   i   l   e   
   p   l   o   t   (   s   l   _   t   i   m   e   2   ,   s   l   _   d   a   t   a   2   (   :   ,   i   d   x   )   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   h   o   l   d       o   n   
   t   i   t   l   e   (   '   I   C   E   6   G   '   )   
   y   l   a   b   e   l   (   '   R   S   L       (   m   )   '   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )                               
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   2       0   ]   ,   '   Y   L   i   m   '   ,   [   -   6   0       1   0   ]   )           
   
   
   n   e   x   t   t   i   l   e   
   
   s   l   r   a   t   e       =       d   i   f   f   (   s   l   _   d   a   t   a   )   .   /   d   i   f   f   (   s   l   _   t   i   m   e   )   ;   
   p   l   o   t   (   t   ,   s   l   r   a   t   e   (   :   ,   i   d   x   )   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   h   o   l   d       o   n   
   
   l   e   g   e   n   d   (   l   i   s   t   (   i   d   x   )   )   
   y   l   a   b   e   l   (   '   R   S   L   R       (   m   m   /   y   r   )   '   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   2       0   ]   ,   '   Y   L   i   m   '   ,   [   -   5       2   0   ]   )   
   
   
   
   
   n   e   x   t   t   i   l   e   
   s   l   r   a   t   e   2       =       d   i   f   f   (   s   l   _   d   a   t   a   2   )   .   /   d   i   f   f   (   s   l   _   t   i   m   e   2   )   ;   
   p   l   o   t   (   t   2   ,   s   l   r   a   t   e   2   (   :   ,   i   d   x   )   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   h   o   l   d       o   n   
   y   l   a   b   e   l   (   '   R   S   L   R       (   m   m   /   y   r   )   '   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   2       0   ]   ,   '   Y   L   i   m   '   ,   [   -   5       2   0   ]   )   
   
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               s   a   v   e   a   s   (   g   c   f   ,   '   F   i   g   S   1   .   s   v   g   '   )   
   
   
   %   %       j   u   s   t       P   o       r   i   v   e   r   
   
   %   %       a   l   t   e   r   n   a   t   i   v   e   
   c   l   r   
   
   l   i   s   t       =       [   "   D   v   i   n   a   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   a   m   p   u   t   r   a   "   ,   "   M   a   c   K   e   n   z   i   e   "   ,   "   M   e   k   o   n   g   "   ,   "   M   i   s   s   i   s   s   i   p   p   i   "   ,   "   N   i   g   e   r   "   ,   "   N   i   l   e   "   ,   "   P   a   r   a   n   a   "   ,   "   W   a   i   p   a   o   a   "   ,   "   Y   e   l   l   o   w   "   ,   "   G   a   n   g   e   s   -   B   r   a   h   m   a   p   u   t   r   a   "   ,   "   S   o   n   g       H   o   n   g   "   ,   "   P   o   "   ,   "   R   h   i   n   e   -   M   e   u   s   e   "   ,   "   T   a   m   a   "   ,   "   T   r   i   n   i   t   y   "   ]   ;   
   
   f       =       '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   '   ;   
   
   t   i   l   e   d   l   a   y   o   u   t   (   '   f   l   o   w   '   )   
   n   e   x   t   t   i   l   e   
   
   %   p   a   s   t       S   L   
   %   d       =       r   e   a   d   t   a   b   l   e   (   [   f       '   a   n   u   _   r   s   l   \   d   e   l   t   a   _   r   s   l   _   a   n   u   h   r   _   7   1   p   2   3   0   .   d   a   t   '   ]   )   ;   
   d       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   i   c   e   6   g   \   d   e   l   t   a   _   r   s   l   _   i   c   e   6   g   _   v   m   5   a   2   .   d   a   t   '   )   ;   
   s   l   _   t   i   m   e       =       d   {   :   ,   1   }   ;   
   s   l   _   d   a   t   a       =       d   {   :   ,   [   2   :   2   :   e   n   d   ]   }   ;   
   
   f   i   d       =       f   o   p   e   n   (   [   f       '   i   c   e   6   g   \   g   m   s   l   _   d   e   f   .   p   l   o   t   '   ]   )   ;       
   g   m   s   l   _   d   e   f       =       t   e   x   t   s   c   a   n   (   f   i   d   ,   '   %   f       %   f   '   )   ;       
   g   m   s   l   _   d   e   f       =       [   g   m   s   l   _   d   e   f   {   :   }   ]   ;   
   f   c   l   o   s   e   (   f   i   d   )   ;   
   
   f   i   d       =       f   o   p   e   n   (   [   f       '   i   c   e   6   g   \   g   m   s   l   _   m   w   .   p   l   o   t   '   ]   )   ;       
   g   m   s   l   _   m   w       =       t   e   x   t   s   c   a   n   (   f   i   d   ,   '   %   f       %   f   '   )   ;       
   g   m   s   l   _   m   w       =       [   g   m   s   l   _   m   w   {   :   }   ]   ;   
   f   c   l   o   s   e   (   f   i   d   )   ;   
   
   f   i   d       =       f   o   p   e   n   (   [   f       '   i   c   e   6   g   \   g   m   s   l   _   n   e   t   .   p   l   o   t   '   ]   )   ;       
   g   m   s   l   _   n   e   t       =       t   e   x   t   s   c   a   n   (   f   i   d   ,   '   %   f       %   f   '   )   ;       
   g   m   s   l   _   n   e   t       =       [   g   m   s   l   _   n   e   t   {   :   }   ]   ;   
   f   c   l   o   s   e   (   f   i   d   )   ;   
   
   %   f   u   t   u   r   e       S   L   
   l   o   a   d   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   d   e   l   t   a   -   t   s   .   m   a   t   '   )   
   
   t       =       s   l   _   t   i   m   e   (   2   :   e   n   d   )   ;   
   
   s   l   r   a   t   e       =       d   i   f   f   (   s   l   _   d   a   t   a   )   .   /   d   i   f   f   (   s   l   _   t   i   m   e   )   ;   
   p   l   o   t   (   t   ,   s   l   r   a   t   e   (   :   ,   [   1   1   :   1   6   ]   )   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   ;       %   ,   '   C   o   l   o   r   '   ,   [   0   .   5       0   .   5       0   .   5   ]   )   
   
   h   o   l   d       o   n   
   
   g   m   s   l       =       d   i   f   f   (   g   m   s   l   _   n   e   t   (   :   ,   2   )   )   .   /   d   i   f   f   (   g   m   s   l   _   n   e   t   (   :   ,   1   )   )   ;   
   p   l   o   t   (   g   m   s   l   _   n   e   t   (   2   :   e   n   d   ,   1   )   ,   g   m   s   l   ,   '   k   '   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   
   m   i       =       [   -   7   ,   -   8   ,   -   7   .   5   ,   -   5   ,   -   7   .   5   ,   -   2   .   5   ]   ;       %   t   i   m   e       i   n       k   a       B   P       o   f       m   i   n   i   m   u   m       i   n       e   x   t   e   n   t       f   r   o   m       q   u   o   c   k   
   [   ~   ,   i   s   a   ]       =       i   s   m   e   m   b   e   r   (   m   i   ,   t   )   ;   
   a       =       f   l   i   p   u   d   (   g   e   t   (   g   c   a   ,   '   C   h   i   l   d   r   e   n   '   )   )   ;       a       =       g   e   t   (   a   ,   '   C   o   l   o   r   '   )   ;   
   f   o   r       i   i   =   1   :   6   ,   
   s   c   a   t   t   e   r   (   t   (   i   s   a   (   i   i   )   )   ,   s   l   r   a   t   e   (   i   s   a   (   i   i   )   ,   1   0   +   i   i   )   ,   5   0   ,   '   M   a   r   k   e   r   F   a   c   e   C   o   l   o   r   '   ,   c   e   l   l   2   m   a   t   (   a   (   i   i   )   )   ,   '   M   a   r   k   e   r   E   d   g   e   C   o   l   o   r   '   ,   '   k   '   )   
   e   n   d   
   
   l   e   g   e   n   d   (   l   i   s   t   (   1   1   :   1   6   )   ,   '   L   o   c   a   t   i   o   n   '   ,   '   S   o   u   t   h   W   e   s   t   '   )   
   
   
   
   s   l   _   f   u   t   u   r   e       =       d   i   f   f   (   [   S   S   P   1   2   6   _   5   0   (   :   ,   4   :   e   n   d   )   ;       S   S   P   2   4   5   _   5   0   (   :   ,   4   :   e   n   d   )   ;   S   S   P   5   8   5   _   5   0   (   :   ,   4   :   e   n   d   )   ]   '   )   .   /   1   0   ;   
   t       =       [   0   ;       d   o   u   b   l   e   (   t   i   m   e   (   1   :   e   n   d   -   1   )   )   .   /   1   0   0   0   -   2   ]   ;   
   %   p   l   o   t   (   t   ,   [   r   e   p   m   a   t   (   s   l   r   a   t   e   (   e   n   d   ,   :   )   ,   1   ,   3   )   ;   s   l   _   f   u   t   u   r   e   (   :   ,   [   1   :   1   0   ,   1   2   :   2   1   ,   2   3   :   3   2   ]   )   ]   ,   '   C   o   l   o   r   '   ,   [   0   .   5       0   .   5       0   .   5   ]   )   
   p   l   o   t   (   t   ,   [   r   e   p   m   a   t   (   g   m   s   l   (   e   n   d   )   ,   1   ,   3   )   ;   s   l   _   f   u   t   u   r   e   (   :   ,   [   1   1       2   2       3   3   ]   )   ]   ,   '   k   '   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   )   
   
   y   l   a   b   e   l   (   '   R   S   L   R       (   m   m   /   y   r   )   '   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   1       0   .   5   ]   )   ;   
   
   n   e   x   t   t   i   l   e   
   %   q   u   o   c   k       d   a   t   a   
   t       =       r   e   a   d   t   a   b   l   e   (   '   D   :   \   D   r   i   v   e   \   g   i   t   h   u   b   \   G   l   o   b   a   l   D   e   l   t   a   C   h   a   n   g   e   \   p   u   b   _   f   i   g   u   r   e   s   \   a   n   n   u   a   l   r   e   v   i   e   w   s   2   0   2   3   \   q   u   o   c   k   .   x   l   s   x   '   )   ;   
   [   l   ,   i   d   x   ]       =       u   n   i   q   u   e   (   t   .   B   a   s   i   n   I   D   2   )   ;   
   f   o   r       i   i   =   1   :   l   e   n   g   t   h   (   l   )   ,   
                   h   o   l   d       o   n   
                   d       =       t   .   n   o   n   d   i   m   e   s   i   o   n   a   l   i   z   e   d   D   i   s   t   a   n   c   e   _   d   i   s   t   a   n   c   e   A   t   T   _   d   i   s   t   a   n   c   e   T   0   _   (   t   .   B   a   s   i   n   I   D   2   =   =   l   (   i   i   )   )   ;   
                   t   t       =       -   t   .   T   i   m   e   _   k   y   r   B   P   (   t   .   B   a   s   i   n   I   D   2   =   =   l   (   i   i   )   )   ;   
                   p   l   o   t   (   t   t   ,   d   *   1   0   0   ,   '   L   i   n   e   W   i   d   t   h   '   ,   2   ,   '   C   o   l   o   r   '   ,   a   {   i   i   }   )   ;   
                   %   [   ~   ,   j   j   ]       =       m   i   n   (   d   )   ;   
                   %   s   c   a   t   t   e   r   (   t   t   (   j   j   )   ,   d   (   j   j   )   .   *   1   0   0   ,   '   M   a   r   k   e   r   F   a   c   e   C   o   l   o   r   '   ,   g   e   t   (   a   ,   '   C   o   l   o   r   '   )   )   
                   %   t   t   (   j   j   )   
                   
   e   n   d                   
   l   e   g   e   n   d   (   t   .   D   e   l   t   a   (   i   d   x   )   ,   '   L   o   c   a   t   i   o   n   '   ,   '   S   o   u   t   h   W   e   s   t   '   )   
   s   e   t   (   g   c   a   ,   '   X   L   i   m   '   ,   [   -   1   1       0   .   5   ]   )   
   x   l   a   b   e   l   (   '   T   i   m   e       (   k   a       B   P   )   '   )   
   y   l   a   b   e   l   (   '   D   e   l   t   a       L   e   n   g   t   h       c   o   m   p   a   r   e   d       t   o       m   o   d   e   r   n       (   %   )   '   )   
   s   e   t   (   g   c   a   ,   '   Y   G   r   i   d   '   ,   '   o   n   '   )   
