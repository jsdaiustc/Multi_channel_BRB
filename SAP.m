function i_sap=SAP(i_b)
i_shift=[i_b(:,2:end),i_b(:,1)];
i_sap=sum(i_b.*i_shift,2);
i_sap=i_sap-mean(i_sap);