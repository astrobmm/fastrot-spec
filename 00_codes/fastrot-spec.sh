rm -f *~
echo "   "
echo "     *****************************************"
echo "     *        PROGRAM  star3d running        *"
echo "     *****************************************"
./star3d.exe
echo "   "
echo "     *****************************************"
echo "     *        PROGRAM  spec3d running        *"
echo "     *****************************************"

./spec3d.exe
echo "   "
echo "     *****************************************"
echo "     *       PROGRAM  normspec running       *"
echo "     *****************************************"
./normspec.exe
mv *.params ../11_results/.
mv *.a      ../11_results/.
mv *.n1     ../11_results/.
echo "   "
echo "     **********************************************"
echo "     *       Synthetic spectra (in physical       *"
echo "     *       units, and normalized) and stellar   *"
echo "     *       parameters are in folder 11_resuts   *"
echo "     **********************************************"
echo "   "

exit
