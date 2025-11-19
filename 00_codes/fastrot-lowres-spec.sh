rm -f *~
echo "   "
echo "     *****************************************"
echo "     *        PROGRAM  star3d running        *"
echo "     *****************************************"
./star3d.exe
./grid10.exe

echo "   "
echo "     *****************************************"
echo "     *      PROGRAM  lowrespec3d running     *"
echo "     *****************************************"

./lowrespec3d.exe
./mergefiles.exe

mv in*.params    ../11_results/.
cp low-res-in*.a ../11_results/.
mv low-res-in*.a ../22_for_plotting/.

mv mer10.a       ../22_for_plotting/.
mv par10.a       ../22_for_plotting/.
mv colat-r-t.a   ../22_for_plotting/.
mv 4plots.a      ../22_for_plotting/.

rm -f *.a

echo "   "
echo "     *******************************************************"
echo "     *                                                     *"  
echo "     *     Synthetic spectrum in physical units and        *"
echo "     *     relevant stellar parameters are in folder:      *"
echo "     *                                                     *"
echo "     *     11_results/                                     *"
echo "     *                                                     *"
echo "     *     Additional output files for plotting,           *"
echo "     *     including the spectra, are in folder:           *"
echo "     *                                                     *"
echo "     *     22_for_plotting/                                *"
echo "     *                                                     *"
echo "     *******************************************************"
echo "   "

exit
