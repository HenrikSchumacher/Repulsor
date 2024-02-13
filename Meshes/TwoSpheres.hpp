// Some hard-coded example featuring two meshed sphere is close vicinity to to each other.
constexpr Int dom_dim      = 2;
constexpr Int amb_dim      = 3;
constexpr Int vertex_count           = 162;
constexpr Int simplex_count          = 320;

constexpr Int obstacle_vertex_count  = 162;
constexpr Int obstacle_simplex_count = 320;

Real vertex_coordinates [vertex_count][amb_dim] = {
    {-0.0000000000000000, 0.0000000000000000, 0.0039062500000000},
    {0.0000000000000000, 0.0000000000000000, 2.0039062500000000},
    {-0.8944271909999160, 0.0000000000000000, 0.5566926545000421},
    {0.8944271909999160, 0.0000000000000000, 1.4511198454999580},
    {0.7236067977499788, -0.5257311121191335, 0.5566926545000421},
    {0.7236067977499788, 0.5257311121191335, 0.5566926545000421},
    {-0.7236067977499788, -0.5257311121191335, 1.4511198454999580},
    {-0.7236067977499788, 0.5257311121191335, 1.4511198454999580},
    {-0.2763932022500210, -0.8506508083520400, 0.5566926545000421},
    {-0.2763932022500210, 0.8506508083520400, 0.5566926545000421},
    {0.2763932022500209, -0.8506508083520400, 1.4511198454999580},
    {0.2763932022500209, 0.8506508083520400, 1.4511198454999580},
    {-0.5257311121191337, 0.0000000000000000, 0.1532554416479601},
    {0.4253254041760199, -0.3090169943749475, 0.1532554416479601},
    {0.4253254041760199, 0.3090169943749475, 0.1532554416479601},
    {-0.1624598481164532, -0.4999999999999999, 0.1532554416479601},
    {-0.1624598481164532, 0.4999999999999999, 0.1532554416479601},
    {0.5257311121191337, 0.0000000000000000, 1.8545570583520400},
    {-0.4253254041760199, -0.3090169943749475, 1.8545570583520400},
    {-0.4253254041760199, 0.3090169943749475, 1.8545570583520400},
    {0.1624598481164532, -0.4999999999999999, 1.8545570583520400},
    {0.1624598481164532, 0.4999999999999999, 1.8545570583520400},
    {-0.9510565162951540, -0.3090169943749474, 1.0039062500000000},
    {-0.9510565162951540, 0.3090169943749474, 1.0039062500000000},
    {-0.6881909602355868, -0.5000000000000001, 0.4781751378808664},
    {-0.6881909602355868, 0.5000000000000001, 0.4781751378808664},
    {0.9510565162951540, -0.3090169943749474, 1.0039062500000000},
    {0.9510565162951540, 0.3090169943749474, 1.0039062500000000},
    {0.6881909602355868, -0.5000000000000001, 1.5296373621191330},
    {0.6881909602355868, 0.5000000000000001, 1.5296373621191330},
    {0.8506508083520400, 0.0000000000000000, 0.4781751378808664},
    {0.2628655560595668, -0.8090169943749480, 0.4781751378808665},
    {0.5877852522924730, -0.8090169943749480, 1.0039062500000000},
    {0.2628655560595668, 0.8090169943749480, 0.4781751378808665},
    {0.5877852522924730, 0.8090169943749480, 1.0039062500000000},
    {-0.8506508083520400, 0.0000000000000000, 1.5296373621191330},
    {-0.5877852522924730, -0.8090169943749480, 1.0039062500000000},
    {-0.2628655560595668, -0.8090169943749480, 1.5296373621191330},
    {-0.5877852522924730, 0.8090169943749480, 1.0039062500000000},
    {-0.2628655560595668, 0.8090169943749480, 1.5296373621191330},
    {0.0000000000000000, -1.0000000000000000, 1.0039062500000000},
    {0.0000000000000000, 1.0000000000000000, 1.0039062500000000},
    {-0.2536652945350961, 0.0000000000000000, 0.0366141956811283},
    {0.2052195341620191, -0.1491007191461559, 0.0366141956811283},
    {0.2052195341620191, 0.1491007191461559, 0.0366141956811283},
    {-0.0783868868944712, -0.2412500313255325, 0.0366141956811283},
    {-0.0783868868944712, 0.2412500313255325, 0.0366141956811283},
    {0.2536652945350961, 0.0000000000000000, 1.9711983043188720},
    {-0.2052195341620191, -0.1491007191461559, 1.9711983043188720},
    {-0.2052195341620191, 0.1491007191461559, 1.9711983043188720},
    {0.0783868868944712, -0.2412500313255325, 1.9711983043188720},
    {0.0783868868944712, 0.2412500313255325, 1.9711983043188720},
    {-0.7517297465983703, 0.0000000000000000, 0.3444349556443246},
    {-0.9569492807603890, -0.1491007191461559, 0.7548740239683629},
    {-0.9569492807603890, 0.1491007191461559, 0.7548740239683629},
    {-0.8301166334928410, -0.2412500313255325, 0.5012087294332669},
    {-0.8301166334928410, 0.2412500313255325, 0.5012087294332669},
    {0.7517297465983703, 0.0000000000000000, 1.6633775443556750},
    {0.9569492807603890, -0.1491007191461559, 1.2529384760316370},
    {0.9569492807603890, 0.1491007191461559, 1.2529384760316370},
    {0.8301166334928410, -0.2412500313255325, 1.5066037705667330},
    {0.8301166334928410, 0.2412500313255325, 1.5066037705667330},
    {0.6081621401752544, -0.4418556587600801, 0.3444349556443245},
    {0.8618274347103510, -0.4418556587600799, 0.7548740239683629},
    {0.8133816743372740, -0.2927549396139240, 0.5012087294332668},
    {0.5297752532807833, -0.6831056900856125, 0.5012087294332669},
    {0.6865490270697254, -0.6831056900856125, 0.7548740239683629},
    {0.6081621401752544, 0.4418556587600801, 0.3444349556443245},
    {0.8618274347103510, 0.4418556587600799, 0.7548740239683629},
    {0.8133816743372740, 0.2927549396139240, 0.5012087294332668},
    {0.5297752532807833, 0.6831056900856125, 0.5012087294332669},
    {0.6865490270697254, 0.6831056900856125, 0.7548740239683629},
    {-0.6081621401752544, -0.4418556587600801, 1.6633775443556760},
    {-0.8618274347103510, -0.4418556587600799, 1.2529384760316370},
    {-0.8133816743372740, -0.2927549396139240, 1.5066037705667330},
    {-0.6865490270697254, -0.6831056900856125, 1.2529384760316370},
    {-0.5297752532807833, -0.6831056900856125, 1.5066037705667330},
    {-0.6081621401752544, 0.4418556587600801, 1.6633775443556760},
    {-0.8618274347103510, 0.4418556587600799, 1.2529384760316370},
    {-0.8133816743372740, 0.2927549396139240, 1.5066037705667330},
    {-0.6865490270697254, 0.6831056900856125, 1.2529384760316370},
    {-0.5297752532807833, 0.6831056900856125, 1.5066037705667330},
    {-0.2322972668760692, -0.7149374739952846, 0.3444349556443246},
    {-0.4859625614111653, -0.7149374739952846, 0.5012087294332669},
    {-0.0270777327140500, -0.8640381931414410, 0.5012087294332669},
    {-0.4375168010380884, -0.8640381931414410, 0.7548740239683629},
    {-0.1539103799815980, -0.9561875053208170, 0.7548740239683629},
    {-0.2322972668760692, 0.7149374739952846, 0.3444349556443246},
    {-0.4859625614111653, 0.7149374739952846, 0.5012087294332669},
    {-0.0270777327140500, 0.8640381931414410, 0.5012087294332669},
    {-0.4375168010380884, 0.8640381931414410, 0.7548740239683629},
    {-0.1539103799815980, 0.9561875053208170, 0.7548740239683629},
    {0.2322972668760692, -0.7149374739952846, 1.6633775443556750},
    {0.4859625614111653, -0.7149374739952846, 1.5066037705667330},
    {0.4375168010380884, -0.8640381931414410, 1.2529384760316370},
    {0.0270777327140500, -0.8640381931414410, 1.5066037705667330},
    {0.1539103799815980, -0.9561875053208170, 1.2529384760316370},
    {0.2322972668760692, 0.7149374739952846, 1.6633775443556750},
    {0.4859625614111653, 0.7149374739952846, 1.5066037705667330},
    {0.4375168010380884, 0.8640381931414410, 1.2529384760316370},
    {0.0270777327140500, 0.8640381931414410, 1.5066037705667330},
    {0.1539103799815980, 0.9561875053208170, 1.2529384760316370},
    {-0.3694077237683695, -0.2683904214913772, 0.1142408548811650},
    {-0.3694077237683695, 0.2683904214913772, 0.1142408548811650},
    {-0.6305367639340622, -0.2683904214913772, 0.2756274771531995},
    {-0.6305367639340622, 0.2683904214913772, 0.2756274771531995},
    {0.4566130579911694, 0.0000000000000000, 0.1142408548811650},
    {0.1411011947727847, -0.4342648242279586, 0.1142408548811650},
    {0.6678708892100330, -0.1534877987546830, 0.2756274771531995},
    {0.3523590259916484, -0.5877526229826416, 0.2756274771531995},
    {0.1411011947727847, 0.4342648242279586, 0.1142408548811650},
    {0.6678708892100330, 0.1534877987546830, 0.2756274771531995},
    {0.3523590259916484, 0.5877526229826416, 0.2756274771531995},
    {-0.4501010349043868, -0.5167388967348582, 0.2756274771531995},
    {0.0604078836367675, -0.6826132994714395, 0.2756274771531995},
    {-0.4501010349043868, 0.5167388967348582, 0.2756274771531995},
    {0.0604078836367675, 0.6826132994714395, 0.2756274771531995},
    {0.3694077237683695, -0.2683904214913772, 1.8935716451188350},
    {0.3694077237683695, 0.2683904214913772, 1.8935716451188350},
    {0.6305367639340622, -0.2683904214913772, 1.7321850228468000},
    {0.6305367639340622, 0.2683904214913772, 1.7321850228468000},
    {-0.4566130579911694, 0.0000000000000000, 1.8935716451188350},
    {-0.1411011947727847, -0.4342648242279586, 1.8935716451188350},
    {-0.6678708892100330, -0.1534877987546830, 1.7321850228468000},
    {-0.3523590259916484, -0.5877526229826416, 1.7321850228468000},
    {-0.1411011947727847, 0.4342648242279586, 1.8935716451188350},
    {-0.6678708892100330, 0.1534877987546830, 1.7321850228468000},
    {-0.3523590259916484, 0.5877526229826416, 1.7321850228468000},
    {0.4501010349043868, -0.5167388967348582, 1.7321850228468000},
    {-0.0604078836367675, -0.6826132994714395, 1.7321850228468000},
    {0.4501010349043868, 0.5167388967348582, 1.7321850228468000},
    {-0.0604078836367675, 0.6826132994714395, 1.7321850228468000},
    {-0.9999444877024320, 0.0000000000000000, 1.0144429246899380},
    {-0.8588432929296470, -0.4342648242279585, 0.7322405351443689},
    {-0.9500732787556030, -0.1534877987546830, 1.2755719648556310},
    {-0.8089720839828180, -0.5877526229826416, 0.9933695753100620},
    {-0.8588432929296470, 0.4342648242279585, 0.7322405351443689},
    {-0.9500732787556030, 0.1534877987546830, 1.2755719648556310},
    {-0.8089720839828180, 0.5877526229826416, 0.9933695753100620},
    {-0.6784075638999713, -0.6826132994714395, 0.7322405351443689},
    {-0.6784075638999713, 0.6826132994714395, 0.7322405351443689},
    {0.9999444877024320, 0.0000000000000000, 0.9933695753100620},
    {0.8588432929296470, -0.4342648242279585, 1.2755719648556310},
    {0.9500732787556030, -0.1534877987546830, 0.7322405351443689},
    {0.8089720839828180, -0.5877526229826416, 1.0144429246899380},
    {0.8588432929296470, 0.4342648242279585, 1.2755719648556310},
    {0.9500732787556030, 0.1534877987546830, 0.7322405351443689},
    {0.8089720839828180, 0.5877526229826416, 1.0144429246899380},
    {0.6784075638999713, -0.6826132994714395, 1.2755719648556310},
    {0.6784075638999713, 0.6826132994714395, 1.2755719648556310},
    {0.4395643602144484, -0.8561430444740190, 0.7322405351443689},
    {0.1476132178595674, -0.9510037209628170, 0.7322405351443689},
    {0.3089998401316020, -0.9510037209628170, 0.9933695753100620},
    {0.4395643602144484, 0.8561430444740190, 0.7322405351443689},
    {0.1476132178595674, 0.9510037209628170, 0.7322405351443689},
    {0.3089998401316020, 0.9510037209628170, 0.9933695753100620},
    {-0.4395643602144484, -0.8561430444740190, 1.2755719648556310},
    {-0.3089998401316020, -0.9510037209628170, 1.0144429246899380},
    {-0.1476132178595674, -0.9510037209628170, 1.2755719648556310},
    {-0.4395643602144484, 0.8561430444740190, 1.2755719648556310},
    {-0.3089998401316020, 0.9510037209628170, 1.0144429246899380},
    {-0.1476132178595674, 0.9510037209628170, 1.2755719648556310}
};

Int simplices [simplex_count][dom_dim+1] = {
    {  5,  68,  69}, { 27, 146,  68}, { 30,  69, 146}, {146,  69,  68}, {  3,  58,  59},
    { 26, 141,  58}, { 27,  59, 141}, {141,  59,  58}, {  4,  64,  63}, { 30, 143,  64},
    { 26,  63, 143}, {143,  63,  64}, { 26, 143, 141}, { 30, 146, 143}, { 27, 141, 146},
    {146, 141, 143}, {  3,  59,  61}, { 27, 145,  59}, { 29,  61, 145}, {145,  61,  59},
    {  5,  71,  68}, { 34, 147,  71}, { 27,  68, 147}, {147,  68,  71}, { 11,  98,  99},
    { 29, 149,  98}, { 34,  99, 149}, {149,  99,  98}, { 34, 149, 147}, { 29, 145, 149},
    { 27, 147, 145}, {145, 147, 149}, { 10,  96,  94}, { 40, 152,  96}, { 32,  94, 152},
    {152,  94,  96}, {  8,  84,  86}, { 31, 151,  84}, { 40,  86, 151}, {151,  86,  84},
    {  4,  66,  65}, { 32, 150,  66}, { 31,  65, 150}, {150,  65,  66}, { 31, 150, 151},
    { 32, 152, 150}, { 40, 151, 152}, {152, 151, 150}, {  8,  86,  85}, { 40, 157,  86},
    {36, 85, 157}, {157, 85, 86}, {10, 95, 96}, {37, 158, 95}, {40, 96, 158},
    {158, 96, 95}, {6, 75, 76}, {36, 156, 75}, {37, 76, 156}, {156, 76, 75},
    {37, 156, 158}, {36, 157, 156}, {40, 158, 157}, {157, 158, 156}, {1, 50, 47},
    {20, 117, 50}, {17, 47, 117}, {117, 47, 50}, {10, 93, 92}, {28, 128, 93},
    {20, 92, 128}, {128, 92, 93}, {3, 57, 60}, {17, 119, 57}, {28, 60, 119},
    {119, 60, 57}, {28, 119, 128}, {17, 117, 119}, {20, 128, 117}, {117, 128, 119},
    {1, 47, 51}, {17, 118, 47}, {21, 51, 118}, {118, 51, 47}, {3, 61, 57},
    {29, 120, 61}, {17, 57, 120}, {120, 57, 61}, {11, 97, 98}, {21, 130, 97},
    {29, 98, 130}, {130, 98, 97}, {29, 130, 120}, {21, 118, 130}, {17, 120, 118},
    {118, 120, 130}, {9, 89, 87}, {33, 116, 89}, {16, 87, 116}, {116, 87, 89},
    {5, 67, 70}, {14, 112, 67}, {33, 70, 112}, {112, 70, 67}, {0, 46, 44},
    {16, 110, 46}, {14, 44, 110}, {110, 44, 46}, {14, 110, 112}, {16, 116, 110},
    {33, 112, 116}, {116, 112, 110}, {5, 69, 67}, {30, 111, 69}, {14, 67, 111},
    {111, 67, 69}, {4, 62, 64}, {13, 108, 62}, {30, 64, 108}, {108, 64, 62},
    {0, 44, 43}, {14, 106, 44}, {13, 43, 106}, {106, 43, 44}, {13, 106, 108},
    {14, 111, 106}, {30, 108, 111}, {111, 108, 106}, {9, 87, 88}, {16, 115, 87},
    {25, 88, 115}, {115, 88, 87}, {0, 42, 46}, {12, 103, 42}, {16, 46, 103},
    {103, 46, 42}, {2, 56, 52}, {25, 105, 56}, {12, 52, 105}, {105, 52, 56},
    {12, 105, 103}, {25, 115, 105}, {16, 103, 115}, {115, 103, 105}, {0, 45, 42},
    {15, 102, 45}, {12, 42, 102}, {102, 42, 45}, {8, 83, 82}, {24, 113, 83},
    {15, 82, 113}, {113, 82, 83}, {2, 52, 55}, {12, 104, 52}, {24, 55, 104},
    {104, 55, 52}, {24, 104, 113}, {12, 102, 104}, {15, 113, 102}, {102, 113, 104},
    {2, 54, 56}, {23, 136, 54}, {25, 56, 136}, {136, 56, 54}, {7, 80, 78},
    {38, 138, 80}, {23, 78, 138}, {138, 78, 80}, {9, 88, 90}, {25, 140, 88},
    {38, 90, 140}, {140, 90, 88}, {38, 140, 138}, {25, 136, 140}, {23, 138, 136},
    {136, 138, 140}, {7, 78, 79}, {23, 137, 78}, {35, 79, 137}, {137, 79, 78},
    {2, 53, 54}, {22, 132, 53}, {23, 54, 132}, {132, 54, 53}, {6, 74, 73},
    {35, 134, 74}, {22, 73, 134}, {134, 73, 74}, {22, 134, 132}, {35, 137, 134},
    {23, 132, 137}, {137, 132, 134}, {11, 100, 97}, {39, 131, 100}, {21, 97, 131},
    {131, 97, 100}, {7, 77, 81}, {19, 127, 77}, {39, 81, 127}, {127, 81, 77},
    {1, 51, 49}, {21, 125, 51}, {19, 49, 125}, {125, 49, 51}, {19, 125, 127},
    {21, 131, 125}, {39, 127, 131}, {131, 127, 125}, {1, 49, 48}, {19, 121, 49},
    {18, 48, 121}, {121, 48, 49}, {7, 79, 77}, {35, 126, 79}, {19, 77, 126},
    {126, 77, 79}, {6, 72, 74}, {18, 123, 72}, {35, 74, 123}, {123, 74, 72},
    {35, 123, 126}, {18, 121, 123}, {19, 126, 121}, {121, 126, 123}, {0, 43, 45},
    {13, 107, 43}, {15, 45, 107}, {107, 45, 43}, {4, 65, 62}, {31, 109, 65},
    {13, 62, 109}, {109, 62, 65}, {8, 82, 84}, {15, 114, 82}, {31, 84, 114},
    {114, 84, 82}, {31, 114, 109}, {15, 107, 114}, {13, 109, 107}, {107, 109, 114},
    {2, 55, 53}, {24, 133, 55}, {22, 53, 133}, {133, 53, 55}, {8, 85, 83},
    {36, 139, 85}, {24, 83, 139}, {139, 83, 85}, {6, 73, 75}, {22, 135, 73},
    {36, 75, 135}, {135, 75, 73}, {36, 135, 139}, {22, 133, 135}, {24, 139, 133},
    {133, 139, 135}, {3, 60, 58}, {28, 142, 60}, {26, 58, 142}, {142, 58, 60},
    {10, 94, 93}, {32, 148, 94}, {28, 93, 148}, {148, 93, 94}, {4, 63, 66},
    {26, 144, 63}, {32, 66, 144}, {144, 66, 63}, {32, 144, 148}, {26, 142, 144},
    {28, 148, 142}, {142, 148, 144}, {10, 92, 95}, {20, 129, 92}, {37, 95, 129},
    {129, 95, 92}, {1, 48, 50}, {18, 122, 48}, {20, 50, 122}, {122, 50, 48},
    {6, 76, 72}, {37, 124, 76}, {18, 72, 124}, {124, 72, 76}, {18, 124, 122},
    {37, 129, 124}, {20, 122, 129}, {129, 122, 124}, {9, 91, 89}, {41, 154, 91},
    {33, 89, 154}, {154, 89, 91}, {11, 99, 101}, {34, 155, 99}, {41, 101, 155},
    {155, 101, 99}, {5, 70, 71}, {33, 153, 70}, {34, 71, 153}, {153, 71, 70},
    {34, 153, 155}, {33, 154, 153}, {41, 155, 154}, {154, 155, 153}, {11, 101, 100},
    {41, 161, 101}, {39, 100, 161}, {161, 100, 101}, {9, 90, 91}, {38, 160, 90},
    {41, 91, 160}, {160, 91, 90}, {7, 81, 80}, {39, 159, 81}, {38, 80, 159},
    {159, 80, 81}, {38, 159, 160}, {39, 161, 159}, {41, 160, 161}, {161, 160, 159}};

Real obstacle_vertex_coordinates [obstacle_vertex_count][3] = {
    {-0.0000000000000000, 0.0000000000000000, -2.0039062500000000},
    {0.0000000000000000, 0.0000000000000000, -0.0039062500000000},
    {-0.8944271909999160, 0.0000000000000000, -1.4511198454999580},
    {0.8944271909999160, 0.0000000000000000, -0.5566926545000421},
    {0.7236067977499788, -0.5257311121191335, -1.4511198454999580},
    {0.7236067977499788, 0.5257311121191335, -1.4511198454999580},
    {-0.7236067977499788, -0.5257311121191335, -0.5566926545000421},
    {-0.7236067977499788, 0.5257311121191335, -0.5566926545000421},
    {-0.2763932022500210, -0.8506508083520400, -1.4511198454999580},
    {-0.2763932022500210, 0.8506508083520400, -1.4511198454999580},
    {0.2763932022500209, -0.8506508083520400, -0.5566926545000421},
    {0.2763932022500209, 0.8506508083520400, -0.5566926545000421},
    {-0.5257311121191337, 0.0000000000000000, -1.8545570583520400},
    {0.4253254041760199, -0.3090169943749475, -1.8545570583520400},
    {0.4253254041760199, 0.3090169943749475, -1.8545570583520400},
    {-0.1624598481164532, -0.4999999999999999, -1.8545570583520400},
    {-0.1624598481164532, 0.4999999999999999, -1.8545570583520400},
    {0.5257311121191337, 0.0000000000000000, -0.1532554416479601},
    {-0.4253254041760199, -0.3090169943749475, -0.1532554416479601},
    {-0.4253254041760199, 0.3090169943749475, -0.1532554416479601},
    {0.1624598481164532, -0.4999999999999999, -0.1532554416479601},
    {0.1624598481164532, 0.4999999999999999, -0.1532554416479601},
    {-0.9510565162951540, -0.3090169943749474, -1.0039062500000000},
    {-0.9510565162951540, 0.3090169943749474, -1.0039062500000000},
    {-0.6881909602355868, -0.5000000000000001, -1.5296373621191330},
    {-0.6881909602355868, 0.5000000000000001, -1.5296373621191330},
    {0.9510565162951540, -0.3090169943749474, -1.0039062500000000},
    {0.9510565162951540, 0.3090169943749474, -1.0039062500000000},
    {0.6881909602355868, -0.5000000000000001, -0.4781751378808664},
    {0.6881909602355868, 0.5000000000000001, -0.4781751378808664},
    {0.8506508083520400, 0.0000000000000000, -1.5296373621191330},
    {0.2628655560595668, -0.8090169943749480, -1.5296373621191330},
    {0.5877852522924730, -0.8090169943749480, -1.0039062500000000},
    {0.2628655560595668, 0.8090169943749480, -1.5296373621191330},
    {0.5877852522924730, 0.8090169943749480, -1.0039062500000000},
    {-0.8506508083520400, 0.0000000000000000, -0.4781751378808664},
    {-0.5877852522924730, -0.8090169943749480, -1.0039062500000000},
    {-0.2628655560595668, -0.8090169943749480, -0.4781751378808665},
    {-0.5877852522924730, 0.8090169943749480, -1.0039062500000000},
    {-0.2628655560595668, 0.8090169943749480, -0.4781751378808665},
    {0.0000000000000000, -1.0000000000000000, -1.0039062500000000},
    {0.0000000000000000, 1.0000000000000000, -1.0039062500000000},
    {-0.2536652945350961, 0.0000000000000000, -1.9711983043188720},
    {0.2052195341620191, -0.1491007191461559, -1.9711983043188720},
    {0.2052195341620191, 0.1491007191461559, -1.9711983043188720},
    {-0.0783868868944712, -0.2412500313255325, -1.9711983043188720},
    {-0.0783868868944712, 0.2412500313255325, -1.9711983043188720},
    {0.2536652945350961, 0.0000000000000000, -0.0366141956811283},
    {-0.2052195341620191, -0.1491007191461559, -0.0366141956811283},
    {-0.2052195341620191, 0.1491007191461559, -0.0366141956811283},
    {0.0783868868944712, -0.2412500313255325, -0.0366141956811283},
    {0.0783868868944712, 0.2412500313255325, -0.0366141956811283},
    {-0.7517297465983703, 0.0000000000000000, -1.6633775443556750},
    {-0.9569492807603890, -0.1491007191461559, -1.2529384760316370},
    {-0.9569492807603890, 0.1491007191461559, -1.2529384760316370},
    {-0.8301166334928410, -0.2412500313255325, -1.5066037705667330},
    {-0.8301166334928410, 0.2412500313255325, -1.5066037705667330},
    {0.7517297465983703, 0.0000000000000000, -0.3444349556443246},
    {0.9569492807603890, -0.1491007191461559, -0.7548740239683629},
    {0.9569492807603890, 0.1491007191461559, -0.7548740239683629},
    {0.8301166334928410, -0.2412500313255325, -0.5012087294332669},
    {0.8301166334928410, 0.2412500313255325, -0.5012087294332669},
    {0.6081621401752544, -0.4418556587600801, -1.6633775443556760},
    {0.8618274347103510, -0.4418556587600799, -1.2529384760316370},
    {0.8133816743372740, -0.2927549396139240, -1.5066037705667330},
    {0.5297752532807833, -0.6831056900856125, -1.5066037705667330},
    {0.6865490270697254, -0.6831056900856125, -1.2529384760316370},
    {0.6081621401752544, 0.4418556587600801, -1.6633775443556760},
    {0.8618274347103510, 0.4418556587600799, -1.2529384760316370},
    {0.8133816743372740, 0.2927549396139240, -1.5066037705667330},
    {0.5297752532807833, 0.6831056900856125, -1.5066037705667330},
    {0.6865490270697254, 0.6831056900856125, -1.2529384760316370},
    {-0.6081621401752544, -0.4418556587600801, -0.3444349556443245},
    {-0.8618274347103510, -0.4418556587600799, -0.7548740239683629},
    {-0.8133816743372740, -0.2927549396139240, -0.5012087294332668},
    {-0.6865490270697254, -0.6831056900856125, -0.7548740239683629},
    {-0.5297752532807833, -0.6831056900856125, -0.5012087294332669},
    {-0.6081621401752544, 0.4418556587600801, -0.3444349556443245},
    {-0.8618274347103510, 0.4418556587600799, -0.7548740239683629},
    {-0.8133816743372740, 0.2927549396139240, -0.5012087294332668},
    {-0.6865490270697254, 0.6831056900856125, -0.7548740239683629},
    {-0.5297752532807833, 0.6831056900856125, -0.5012087294332669},
    {-0.2322972668760692, -0.7149374739952846, -1.6633775443556750},
    {-0.4859625614111653, -0.7149374739952846, -1.5066037705667330},
    {-0.0270777327140500, -0.8640381931414410, -1.5066037705667330},
    {-0.4375168010380884, -0.8640381931414410, -1.2529384760316370},
    {-0.1539103799815980, -0.9561875053208170, -1.2529384760316370},
    {-0.2322972668760692, 0.7149374739952846, -1.6633775443556750},
    {-0.4859625614111653, 0.7149374739952846, -1.5066037705667330},
    {-0.0270777327140500, 0.8640381931414410, -1.5066037705667330},
    {-0.4375168010380884, 0.8640381931414410, -1.2529384760316370},
    {-0.1539103799815980, 0.9561875053208170, -1.2529384760316370},
    {0.2322972668760692, -0.7149374739952846, -0.3444349556443246},
    {0.4859625614111653, -0.7149374739952846, -0.5012087294332669},
    {0.4375168010380884, -0.8640381931414410, -0.7548740239683629},
    {0.0270777327140500, -0.8640381931414410, -0.5012087294332669},
    {0.1539103799815980, -0.9561875053208170, -0.7548740239683629},
    {0.2322972668760692, 0.7149374739952846, -0.3444349556443246},
    {0.4859625614111653, 0.7149374739952846, -0.5012087294332669},
    {0.4375168010380884, 0.8640381931414410, -0.7548740239683629},
    {0.0270777327140500, 0.8640381931414410, -0.5012087294332669},
    {0.1539103799815980, 0.9561875053208170, -0.7548740239683629},
    {-0.3694077237683695, -0.2683904214913772, -1.8935716451188350},
    {-0.3694077237683695, 0.2683904214913772, -1.8935716451188350},
    {-0.6305367639340622, -0.2683904214913772, -1.7321850228468000},
    {-0.6305367639340622, 0.2683904214913772, -1.7321850228468000},
    {0.4566130579911694, 0.0000000000000000, -1.8935716451188350},
    {0.1411011947727847, -0.4342648242279586, -1.8935716451188350},
    {0.6678708892100330, -0.1534877987546830, -1.7321850228468000},
    {0.3523590259916484, -0.5877526229826416, -1.7321850228468000},
    {0.1411011947727847, 0.4342648242279586, -1.8935716451188350},
    {0.6678708892100330, 0.1534877987546830, -1.7321850228468000},
    {0.3523590259916484, 0.5877526229826416, -1.7321850228468000},
    {-0.4501010349043868, -0.5167388967348582, -1.7321850228468000},
    {0.0604078836367675, -0.6826132994714395, -1.7321850228468000},
    {-0.4501010349043868, 0.5167388967348582, -1.7321850228468000},
    {0.0604078836367675, 0.6826132994714395, -1.7321850228468000},
    {0.3694077237683695, -0.2683904214913772, -0.1142408548811650},
    {0.3694077237683695, 0.2683904214913772, -0.1142408548811650},
    {0.6305367639340622, -0.2683904214913772, -0.2756274771531995},
    {0.6305367639340622, 0.2683904214913772, -0.2756274771531995},
    {-0.4566130579911694, 0.0000000000000000, -0.1142408548811650},
    {-0.1411011947727847, -0.4342648242279586, -0.1142408548811650},
    {-0.6678708892100330, -0.1534877987546830, -0.2756274771531995},
    {-0.3523590259916484, -0.5877526229826416, -0.2756274771531995},
    {-0.1411011947727847, 0.4342648242279586, -0.1142408548811650},
    {-0.6678708892100330, 0.1534877987546830, -0.2756274771531995},
    {-0.3523590259916484, 0.5877526229826416, -0.2756274771531995},
    {0.4501010349043868, -0.5167388967348582, -0.2756274771531995},
    {-0.0604078836367675, -0.6826132994714395, -0.2756274771531995},
    {0.4501010349043868, 0.5167388967348582, -0.2756274771531995},
    {-0.0604078836367675, 0.6826132994714395, -0.2756274771531995},
    {-0.9999444877024320, 0.0000000000000000, -0.9933695753100620},
    {-0.8588432929296470, -0.4342648242279585, -1.2755719648556310},
    {-0.9500732787556030, -0.1534877987546830, -0.7322405351443689},
    {-0.8089720839828180, -0.5877526229826416, -1.0144429246899380},
    {-0.8588432929296470, 0.4342648242279585, -1.2755719648556310},
    {-0.9500732787556030, 0.1534877987546830, -0.7322405351443689},
    {-0.8089720839828180, 0.5877526229826416, -1.0144429246899380},
    {-0.6784075638999713, -0.6826132994714395, -1.2755719648556310},
    {-0.6784075638999713, 0.6826132994714395, -1.2755719648556310},
    {0.9999444877024320, 0.0000000000000000, -1.0144429246899380},
    {0.8588432929296470, -0.4342648242279585, -0.7322405351443689},
    {0.9500732787556030, -0.1534877987546830, -1.2755719648556310},
    {0.8089720839828180, -0.5877526229826416, -0.9933695753100620},
    {0.8588432929296470, 0.4342648242279585, -0.7322405351443689},
    {0.9500732787556030, 0.1534877987546830, -1.2755719648556310},
    {0.8089720839828180, 0.5877526229826416, -0.9933695753100620},
    {0.6784075638999713, -0.6826132994714395, -0.7322405351443689},
    {0.6784075638999713, 0.6826132994714395, -0.7322405351443689},
    {0.4395643602144484, -0.8561430444740190, -1.2755719648556310},
    {0.1476132178595674, -0.9510037209628170, -1.2755719648556310},
    {0.3089998401316020, -0.9510037209628170, -1.0144429246899380},
    {0.4395643602144484, 0.8561430444740190, -1.2755719648556310},
    {0.1476132178595674, 0.9510037209628170, -1.2755719648556310},
    {0.3089998401316020, 0.9510037209628170, -1.0144429246899380},
    {-0.4395643602144484, -0.8561430444740190, -0.7322405351443689},
    {-0.3089998401316020, -0.9510037209628170, -0.9933695753100620},
    {-0.1476132178595674, -0.9510037209628170, -0.7322405351443689},
    {-0.4395643602144484, 0.8561430444740190, -0.7322405351443689},
    {-0.3089998401316020, 0.9510037209628170, -0.9933695753100620},
    {-0.1476132178595674, 0.9510037209628170, -0.7322405351443689}
};

Int obstacle_simplices [obstacle_simplex_count][dom_dim+1] = {
    {  5,  68,  69}, { 27, 146,  68}, { 30,  69, 146}, {146,  69,  68}, {  3,  58,  59},
    { 26, 141,  58}, { 27,  59, 141}, {141,  59,  58}, {  4,  64,  63}, { 30, 143,  64},
    { 26,  63, 143}, {143,  63,  64}, { 26, 143, 141}, { 30, 146, 143}, { 27, 141, 146},
    {146, 141, 143}, {  3,  59,  61}, { 27, 145,  59}, { 29,  61, 145}, {145,  61,  59},
    {  5,  71,  68}, { 34, 147,  71}, { 27,  68, 147}, {147,  68,  71}, { 11,  98,  99},
    { 29, 149,  98}, { 34,  99, 149}, {149,  99,  98}, { 34, 149, 147}, { 29, 145, 149},
    { 27, 147, 145}, {145, 147, 149}, { 10,  96,  94}, { 40, 152,  96}, { 32,  94, 152},
    {152,  94,  96}, {  8,  84,  86}, { 31, 151,  84}, { 40,  86, 151}, {151,  86,  84},
    {  4,  66,  65}, { 32, 150,  66}, { 31,  65, 150}, {150,  65,  66}, { 31, 150, 151},
    { 32, 152, 150}, { 40, 151, 152}, {152, 151, 150}, {  8,  86,  85}, { 40, 157,  86},
    {36, 85, 157}, {157, 85, 86}, {10, 95, 96}, {37, 158, 95}, {40, 96, 158},
    {158, 96, 95}, {6, 75, 76}, {36, 156, 75}, {37, 76, 156}, {156, 76, 75},
    {37, 156, 158}, {36, 157, 156}, {40, 158, 157}, {157, 158, 156}, {1, 50, 47},
    {20, 117, 50}, {17, 47, 117}, {117, 47, 50}, {10, 93, 92}, {28, 128, 93},
    {20, 92, 128}, {128, 92, 93}, {3, 57, 60}, {17, 119, 57}, {28, 60, 119},
    {119, 60, 57}, {28, 119, 128}, {17, 117, 119}, {20, 128, 117}, {117, 128, 119},
    {1, 47, 51}, {17, 118, 47}, {21, 51, 118}, {118, 51, 47}, {3, 61, 57},
    {29, 120, 61}, {17, 57, 120}, {120, 57, 61}, {11, 97, 98}, {21, 130, 97},
    {29, 98, 130}, {130, 98, 97}, {29, 130, 120}, {21, 118, 130}, {17, 120, 118},
    {118, 120, 130}, {9, 89, 87}, {33, 116, 89}, {16, 87, 116}, {116, 87, 89},
    {5, 67, 70}, {14, 112, 67}, {33, 70, 112}, {112, 70, 67}, {0, 46, 44},
    {16, 110, 46}, {14, 44, 110}, {110, 44, 46}, {14, 110, 112}, {16, 116, 110},
    {33, 112, 116}, {116, 112, 110}, {5, 69, 67}, {30, 111, 69}, {14, 67, 111},
    {111, 67, 69}, {4, 62, 64}, {13, 108, 62}, {30, 64, 108}, {108, 64, 62},
    {0, 44, 43}, {14, 106, 44}, {13, 43, 106}, {106, 43, 44}, {13, 106, 108},
    {14, 111, 106}, {30, 108, 111}, {111, 108, 106}, {9, 87, 88}, {16, 115, 87},
    {25, 88, 115}, {115, 88, 87}, {0, 42, 46}, {12, 103, 42}, {16, 46, 103},
    {103, 46, 42}, {2, 56, 52}, {25, 105, 56}, {12, 52, 105}, {105, 52, 56},
    {12, 105, 103}, {25, 115, 105}, {16, 103, 115}, {115, 103, 105}, {0, 45, 42},
    {15, 102, 45}, {12, 42, 102}, {102, 42, 45}, {8, 83, 82}, {24, 113, 83},
    {15, 82, 113}, {113, 82, 83}, {2, 52, 55}, {12, 104, 52}, {24, 55, 104},
    {104, 55, 52}, {24, 104, 113}, {12, 102, 104}, {15, 113, 102}, {102, 113, 104},
    {2, 54, 56}, {23, 136, 54}, {25, 56, 136}, {136, 56, 54}, {7, 80, 78},
    {38, 138, 80}, {23, 78, 138}, {138, 78, 80}, {9, 88, 90}, {25, 140, 88},
    {38, 90, 140}, {140, 90, 88}, {38, 140, 138}, {25, 136, 140}, {23, 138, 136},
    {136, 138, 140}, {7, 78, 79}, {23, 137, 78}, {35, 79, 137}, {137, 79, 78},
    {2, 53, 54}, {22, 132, 53}, {23, 54, 132}, {132, 54, 53}, {6, 74, 73},
    {35, 134, 74}, {22, 73, 134}, {134, 73, 74}, {22, 134, 132}, {35, 137, 134},
    {23, 132, 137}, {137, 132, 134}, {11, 100, 97}, {39, 131, 100}, {21, 97, 131},
    {131, 97, 100}, {7, 77, 81}, {19, 127, 77}, {39, 81, 127}, {127, 81, 77},
    {1, 51, 49}, {21, 125, 51}, {19, 49, 125}, {125, 49, 51}, {19, 125, 127},
    {21, 131, 125}, {39, 127, 131}, {131, 127, 125}, {1, 49, 48}, {19, 121, 49},
    {18, 48, 121}, {121, 48, 49}, {7, 79, 77}, {35, 126, 79}, {19, 77, 126},
    {126, 77, 79}, {6, 72, 74}, {18, 123, 72}, {35, 74, 123}, {123, 74, 72},
    {35, 123, 126}, {18, 121, 123}, {19, 126, 121}, {121, 126, 123}, {0, 43, 45},
    {13, 107, 43}, {15, 45, 107}, {107, 45, 43}, {4, 65, 62}, {31, 109, 65},
    {13, 62, 109}, {109, 62, 65}, {8, 82, 84}, {15, 114, 82}, {31, 84, 114},
    {114, 84, 82}, {31, 114, 109}, {15, 107, 114}, {13, 109, 107}, {107, 109, 114},
    {2, 55, 53}, {24, 133, 55}, {22, 53, 133}, {133, 53, 55}, {8, 85, 83},
    {36, 139, 85}, {24, 83, 139}, {139, 83, 85}, {6, 73, 75}, {22, 135, 73},
    {36, 75, 135}, {135, 75, 73}, {36, 135, 139}, {22, 133, 135}, {24, 139, 133},
    {133, 139, 135}, {3, 60, 58}, {28, 142, 60}, {26, 58, 142}, {142, 58, 60},
    {10, 94, 93}, {32, 148, 94}, {28, 93, 148}, {148, 93, 94}, {4, 63, 66},
    {26, 144, 63}, {32, 66, 144}, {144, 66, 63}, {32, 144, 148}, {26, 142, 144},
    {28, 148, 142}, {142, 148, 144}, {10, 92, 95}, {20, 129, 92}, {37, 95, 129},
    {129, 95, 92}, {1, 48, 50}, {18, 122, 48}, {20, 50, 122}, {122, 50, 48},
    {6, 76, 72}, {37, 124, 76}, {18, 72, 124}, {124, 72, 76}, {18, 124, 122},
    {37, 129, 124}, {20, 122, 129}, {129, 122, 124}, {9, 91, 89}, {41, 154, 91},
    {33, 89, 154}, {154, 89, 91}, {11, 99, 101}, {34, 155, 99}, {41, 101, 155},
    {155, 101, 99}, {5, 70, 71}, {33, 153, 70}, {34, 71, 153}, {153, 71, 70},
    {34, 153, 155}, {33, 154, 153}, {41, 155, 154}, {154, 155, 153}, {11, 101, 100},
    {41, 161, 101}, {39, 100, 161}, {161, 100, 101}, {9, 90, 91}, {38, 160, 90},
    {41, 91, 160}, {160, 91, 90}, {7, 81, 80}, {39, 159, 81}, {38, 80, 159},
    {159, 80, 81}, {38, 159, 160}, {39, 161, 159}, {41, 160, 161}, {161, 160, 159}};