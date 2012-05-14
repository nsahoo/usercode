float cutrefbdt_hzz[2][3] = 
  { 
    { 0.470, 0.004, 0.295 }, // 5 < pt < 10 GeV
    { 0.500, 0.120, 0.600 } // pt > 10 GeV
  };

// special cuts for HZZ
float cutbdt_hzz[2][3][2] = 
  { 
    {// 5 < pt < 10 GeV
      {  0.369, 0.093 }, // |eta| < 0.8
      { -0.025, 0.451 }, // 0.8 < |eta| < 1.479
      {  0.531, 0.595 } // 1.479 < |eta| < 2.5
    },
    { // pt > 10 GeV
      { 0.735, 0.881 },
      { 0.467, 0.731 },
      { 0.795, 0.891 }
    }
  };

float cutmvaiso_hzz[2][3][2] = 
  {
    { // 5 < pt < 10 GeV
      {  0.385,  0.553 }, // |eta| < 0.8
      { -0.083, -0.237 }, // 0.8 < |eta| < 1.479
      { -0.573, -0.573 }  // 1.479 < |eta| < 2.5
    },
    { // pt > 10 GeV
      { 0.413, 0.521 }, // |eta| < 0.8
      { 0.271, 0.531 }, // 0.8 < |eta| < 1.479
      { 0.135, 0.493 }  // 1.479 < |eta| < 2.5
    }
  };
