rm -f *~
echo "   "
echo "     *****************************************"
echo "     *        PROGRAM  star3d running        *"
echo "     *****************************************"
./star3d.exe
./grid10.exe

mv in*.params    ../11_results/.
mv mer*.a        ../22_for_plotting/.
mv par*.a        ../22_for_plotting/.
mv colat-r-t.a   ../22_for_plotting/.

echo "   "
echo "     ******************************************"
echo "     *      PROGRAM  highrespec3d running     *"
echo "     ******************************************"

./highrespec3d.exe
./mergefiles.exe
./normspec.exe

cp in*.a      ../11_results/.
cp in*.n1     ../11_results/.
mv in*.a      ../22_for_plotting/.
mv in*.n1     ../22_for_plotting/.
mv 4plots.a   ../22_for_plotting/.

rm -f *.a

echo "   "
echo "     **************************************************************"
echo "     *                                                            *"  
echo "     *     Synthetic spectrum in physical units (extension .a)    *"
echo "     *     and normalized (extension .n1), and relevant stellar   *"
echo "     *     parameters are in folder:                              *"
echo "     *                                                            *"
echo "     *     11_results/                                            *"
echo "     *                                                            *"
echo "     *     Additional output files for plotting, including the    *"
echo "     *     the spectra, are in folder:                            *"
echo "     *                                                            *"
echo "     *     22_for_plotting/                                       *"
echo "     *                                                            *"
echo "     **************************************************************"
echo "   "

exit
