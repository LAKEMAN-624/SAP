      function MySum(x,y)
       implicit none
      !DEC$ ATTRIBUTES DLLEXPORT :: MySum
      !DEC$ ATTRIBUTES ALIAS:'MySum'::Mysum
      integer x,y,MySum
      MySum=x+y
      end function
