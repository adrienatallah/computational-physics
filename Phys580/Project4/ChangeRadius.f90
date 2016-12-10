Program ChangeRadius

 Implicit None

 Integer :: Radius
 
 Open (unit = 9, file = 'Radius')
 Read (9, *) Radius
 Rewind (9)
 Write (9, *) Radius + 2
 Close(9)
 
 End Program ChangeRadius
