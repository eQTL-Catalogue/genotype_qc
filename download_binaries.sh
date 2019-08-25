#Make bin directory
mkdir bin
cd bin/

#Download GenotypeHarmonizer
wget http://www.molgenis.org/downloads/GenotypeHarmonizer/GenotypeHarmonizer-1.4.20-dist.tar.gz
tar xzfv GenotypeHarmonizer-1.4.20-dist.tar.gz
mv GenotypeHarmonizer-1.4.20-SNAPSHOT/ GenotypeHarmonizer-1.4.20
rm GenotypeHarmonizer-1.4.20-dist.tar.gz
cd ..