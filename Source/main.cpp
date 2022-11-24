#include "stdafx.h"
#include "RandNum.h"
#include "BCH_BPSK.h"

CLCRandNum rndGen0;
CWHRandNum rndGen1;
CLCRandNum rndGen2;
CLCRandNum rndGen3;
CBCH_BPSK *simulator;

int main()
{
	int flag0, flag1, flag2, flag3;
	int GII_Running;
	FILE *fp;
	char temp_str[80] = { ' ' };

	fprintf(stdout, "\n\nRead SeedSetting.txt file.\n");
	if ((fp = fopen("SeedSetting.txt", "r")) == NULL)
	{
		fprintf(stdout, "Cannot open SeedSetting.txt file./n");
		
		return 1;
	}

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d %d %d %d", &flag0, &flag1, &flag2, &flag3);

	fscanf(fp, "%s", temp_str);
	if (flag0 == 2)
		rndGen0.SetSeed(fp);
	else
		rndGen0.SetSeed(flag0);

	fscanf(fp, "%s", temp_str);
	if (flag1 == 2)
		rndGen1.SetSeed(fp);
	else
		rndGen1.SetSeed(flag1);

	fscanf(fp, "%s", temp_str);
	if (flag2 == 2)
		rndGen2.SetSeed(fp);
	else
		rndGen2.SetSeed(flag2);

	fscanf(fp, "%s", temp_str);
	if (flag3 == 2)
		rndGen3.SetSeed(fp);
	else
		rndGen3.SetSeed(flag3);

	fclose(fp);


	simulator = new CBCH_BPSK;
	simulator->Simulator(); //Here!
		
	delete simulator;


	printf("\n\n**********************************************************");
	printf("\nThe simulation is done!");
	printf("\n**********************************************************\n\n");
	system("pause");

	return 0;
}

