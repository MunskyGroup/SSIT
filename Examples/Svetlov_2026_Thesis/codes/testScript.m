eic1 = MBExperimentInputConfigurable();
eic1.InputName = "IDex";
eic1.Values = ["0" "10" "100"];

epc1 = MBExperimentParameterConfigurable();
epc1.ParameterName = "tptlTime";
epc1.Values = [0 25 50 75 100];

tc1 = MBExperimentTimeConfigurable();
tc1.Values = 0:10:180;

ntc1 = MBExperimentNonTimeConfiguration();
ntc1.NonTimeConfigurables = [eic1 epc1];

ntc2 = MBExperimentNonTimeConfiguration();
ntc2.NonTimeConfigurables = [epc1 eic1];

ec1 = MBExperimentConfiguration();
ec1.NonTimeConfiguration = ntc1;
ec1.TimeConfigurable = tc1;

ec2 = MBExperimentConfiguration();
ec2.NonTimeConfiguration = ntc2;
ec2.TimeConfigurable = tc1;