rm -f *~
echo "   "
echo "     *****************************************"
echo "     *        PROGRAM  star3d running        *"
echo "     *****************************************"
./star3d.exe
echo "   "
echo "     *****************************************"
echo "     *      PROGRAM  highrespec3d running    *"
echo "     *****************************************"

./lowrespec3d.exe

mv *.params ../11_results/.
mv *.a      ../11_results/.

echo "   "
echo "     **********************************************"
echo "     *       Synthetic spectrum (in physical      *"
echo "     *       and stellar parameters are in        *"
echo "     *       parameters are in folder 11_results  *"
echo "     **********************************************"
echo "   "

exit
