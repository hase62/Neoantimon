package utility;

import java.util.*;

import processing.core.PApplet;

public class MyColor {
	private static Map <String, int[]> colName2RGBdec =  new HashMap<String, int[]>();
	private static Map <String, Integer>  colName2RGBhex = new HashMap<String, Integer>();
	
	static{
		//each entry was copied from tables in http://web.njit.edu/~kevin/rgb.html,
		//and has Credits "X" (MIT's Xconsortum red/green/blue (RGB) color specifications).
		{int[] tmp = {190,190,190};
		colName2RGBdec.put("GREY", tmp);
		colName2RGBhex.put("GREY", 0xFFBEBEBE);}
		{int[] tmp = {211,211,211};
		colName2RGBdec.put("LIGHTGRAY", tmp);
		colName2RGBhex.put("LIGHTGRAY", 0xFFD3D3D3);}
		{int[] tmp = {119,136,153};
		colName2RGBdec.put("LIGHTSLATEGREY", tmp);
		colName2RGBhex.put("LIGHTSLATEGREY", 0xFF778899);}
		{int[] tmp = {112,128,144};
		colName2RGBdec.put("SLATEGRAY", tmp);
		colName2RGBhex.put("SLATEGRAY", 0xFF708090);}
		{int[] tmp = {198,226,255};
		colName2RGBdec.put("SLATEGRAY1", tmp);
		colName2RGBhex.put("SLATEGRAY1", 0xFFC6E2FF);}
		{int[] tmp = {185,211,238};
		colName2RGBdec.put("SLATEGRAY2", tmp);
		colName2RGBhex.put("SLATEGRAY2", 0xFFB9D3EE);}
		{int[] tmp = {159,182,205};
		colName2RGBdec.put("SLATEGRAY3", tmp);
		colName2RGBhex.put("SLATEGRAY3", 0xFF9FB6CD);}
		{int[] tmp = {108,123,139};
		colName2RGBdec.put("SLATEGRAY4", tmp);
		colName2RGBhex.put("SLATEGRAY4", 0xFF6C7B8B);}
		{int[] tmp = {0,0,0};
		colName2RGBdec.put("BLACK", tmp);
		colName2RGBhex.put("BLACK", 0xFF000000);}
		{int[] tmp = {0,0,0};
		colName2RGBdec.put("GREY0", tmp);
		colName2RGBhex.put("GREY0", 0xFF000000);}
		{int[] tmp = {3,3,3};
		colName2RGBdec.put("GREY1", tmp);
		colName2RGBhex.put("GREY1", 0xFF030303);}
		{int[] tmp = {5,5,5};
		colName2RGBdec.put("GREY2", tmp);
		colName2RGBhex.put("GREY2", 0xFF050505);}
		{int[] tmp = {8,8,8};
		colName2RGBdec.put("GREY3", tmp);
		colName2RGBhex.put("GREY3", 0xFF080808);}
		{int[] tmp = {10,10,10};
		colName2RGBdec.put("GREY4", tmp);
		colName2RGBhex.put("GREY4", 0xFF0A0A0A);}
		{int[] tmp = {13,13,13};
		colName2RGBdec.put("GREY5", tmp);
		colName2RGBhex.put("GREY5", 0xFF0D0D0D);}
		{int[] tmp = {15,15,15};
		colName2RGBdec.put("GREY6", tmp);
		colName2RGBhex.put("GREY6", 0xFF0F0F0F);}
		{int[] tmp = {18,18,18};
		colName2RGBdec.put("GREY7", tmp);
		colName2RGBhex.put("GREY7", 0xFF121212);}
		{int[] tmp = {20,20,20};
		colName2RGBdec.put("GREY8", tmp);
		colName2RGBhex.put("GREY8", 0xFF141414);}
		{int[] tmp = {23,23,23};
		colName2RGBdec.put("GREY9", tmp);
		colName2RGBhex.put("GREY9", 0xFF171717);}
		{int[] tmp = {26,26,26};
		colName2RGBdec.put("GREY10", tmp);
		colName2RGBhex.put("GREY10", 0xFF1A1A1A);}
		{int[] tmp = {28,28,28};
		colName2RGBdec.put("GREY11", tmp);
		colName2RGBhex.put("GREY11", 0xFF1C1C1C);}
		{int[] tmp = {31,31,31};
		colName2RGBdec.put("GREY12", tmp);
		colName2RGBhex.put("GREY12", 0xFF1F1F1F);}
		{int[] tmp = {33,33,33};
		colName2RGBdec.put("GREY13", tmp);
		colName2RGBhex.put("GREY13", 0xFF212121);}
		{int[] tmp = {36,36,36};
		colName2RGBdec.put("GREY14", tmp);
		colName2RGBhex.put("GREY14", 0xFF242424);}
		{int[] tmp = {38,38,38};
		colName2RGBdec.put("GREY15", tmp);
		colName2RGBhex.put("GREY15", 0xFF262626);}
		{int[] tmp = {41,41,41};
		colName2RGBdec.put("GREY16", tmp);
		colName2RGBhex.put("GREY16", 0xFF292929);}
		{int[] tmp = {43,43,43};
		colName2RGBdec.put("GREY17", tmp);
		colName2RGBhex.put("GREY17", 0xFF2B2B2B);}
		{int[] tmp = {46,46,46};
		colName2RGBdec.put("GREY18", tmp);
		colName2RGBhex.put("GREY18", 0xFF2E2E2E);}
		{int[] tmp = {48,48,48};
		colName2RGBdec.put("GREY19", tmp);
		colName2RGBhex.put("GREY19", 0xFF303030);}
		{int[] tmp = {51,51,51};
		colName2RGBdec.put("GREY20", tmp);
		colName2RGBhex.put("GREY20", 0xFF333333);}
		{int[] tmp = {54,54,54};
		colName2RGBdec.put("GREY21", tmp);
		colName2RGBhex.put("GREY21", 0xFF363636);}
		{int[] tmp = {56,56,56};
		colName2RGBdec.put("GREY22", tmp);
		colName2RGBhex.put("GREY22", 0xFF383838);}
		{int[] tmp = {59,59,59};
		colName2RGBdec.put("GREY23", tmp);
		colName2RGBhex.put("GREY23", 0xFF3B3B3B);}
		{int[] tmp = {61,61,61};
		colName2RGBdec.put("GREY24", tmp);
		colName2RGBhex.put("GREY24", 0xFF3D3D3D);}
		{int[] tmp = {64,64,64};
		colName2RGBdec.put("GREY25", tmp);
		colName2RGBhex.put("GREY25", 0xFF404040);}
		{int[] tmp = {66,66,66};
		colName2RGBdec.put("GREY26", tmp);
		colName2RGBhex.put("GREY26", 0xFF424242);}
		{int[] tmp = {69,69,69};
		colName2RGBdec.put("GREY27", tmp);
		colName2RGBhex.put("GREY27", 0xFF454545);}
		{int[] tmp = {71,71,71};
		colName2RGBdec.put("GREY28", tmp);
		colName2RGBhex.put("GREY28", 0xFF474747);}
		{int[] tmp = {74,74,74};
		colName2RGBdec.put("GREY29", tmp);
		colName2RGBhex.put("GREY29", 0xFF4A4A4A);}
		{int[] tmp = {77,77,77};
		colName2RGBdec.put("GREY30", tmp);
		colName2RGBhex.put("GREY30", 0xFF4D4D4D);}
		{int[] tmp = {79,79,79};
		colName2RGBdec.put("GREY31", tmp);
		colName2RGBhex.put("GREY31", 0xFF4F4F4F);}
		{int[] tmp = {82,82,82};
		colName2RGBdec.put("GREY32", tmp);
		colName2RGBhex.put("GREY32", 0xFF525252);}
		{int[] tmp = {84,84,84};
		colName2RGBdec.put("GREY33", tmp);
		colName2RGBhex.put("GREY33", 0xFF545454);}
		{int[] tmp = {87,87,87};
		colName2RGBdec.put("GREY34", tmp);
		colName2RGBhex.put("GREY34", 0xFF575757);}
		{int[] tmp = {89,89,89};
		colName2RGBdec.put("GREY35", tmp);
		colName2RGBhex.put("GREY35", 0xFF595959);}
		{int[] tmp = {92,92,92};
		colName2RGBdec.put("GREY36", tmp);
		colName2RGBhex.put("GREY36", 0xFF5C5C5C);}
		{int[] tmp = {94,94,94};
		colName2RGBdec.put("GREY37", tmp);
		colName2RGBhex.put("GREY37", 0xFF5E5E5E);}
		{int[] tmp = {97,97,97};
		colName2RGBdec.put("GREY38", tmp);
		colName2RGBhex.put("GREY38", 0xFF616161);}
		{int[] tmp = {99,99,99};
		colName2RGBdec.put("GREY39", tmp);
		colName2RGBhex.put("GREY39", 0xFF636363);}
		{int[] tmp = {102,102,102};
		colName2RGBdec.put("GREY40", tmp);
		colName2RGBhex.put("GREY40", 0xFF666666);}
		{int[] tmp = {105,105,105};
		colName2RGBdec.put("GREY41", tmp);
		colName2RGBhex.put("GREY41", 0xFF696969);}
		{int[] tmp = {105,105,105};
		colName2RGBdec.put("DIMGREY", tmp);
		colName2RGBhex.put("DIMGREY", 0xFF696969);}
		{int[] tmp = {107,107,107};
		colName2RGBdec.put("GREY42", tmp);
		colName2RGBhex.put("GREY42", 0xFF6B6B6B);}
		{int[] tmp = {110,110,110};
		colName2RGBdec.put("GREY43", tmp);
		colName2RGBhex.put("GREY43", 0xFF6E6E6E);}
		{int[] tmp = {112,112,112};
		colName2RGBdec.put("GREY44", tmp);
		colName2RGBhex.put("GREY44", 0xFF707070);}
		{int[] tmp = {115,115,115};
		colName2RGBdec.put("GREY45", tmp);
		colName2RGBhex.put("GREY45", 0xFF737373);}
		{int[] tmp = {117,117,117};
		colName2RGBdec.put("GREY46", tmp);
		colName2RGBhex.put("GREY46", 0xFF757575);}
		{int[] tmp = {120,120,120};
		colName2RGBdec.put("GREY47", tmp);
		colName2RGBhex.put("GREY47", 0xFF787878);}
		{int[] tmp = {122,122,122};
		colName2RGBdec.put("GREY48", tmp);
		colName2RGBhex.put("GREY48", 0xFF7A7A7A);}
		{int[] tmp = {125,125,125};
		colName2RGBdec.put("GREY49", tmp);
		colName2RGBhex.put("GREY49", 0xFF7D7D7D);}
		{int[] tmp = {127,127,127};
		colName2RGBdec.put("GREY50", tmp);
		colName2RGBhex.put("GREY50", 0xFF7F7F7F);}
		{int[] tmp = {130,130,130};
		colName2RGBdec.put("GREY51", tmp);
		colName2RGBhex.put("GREY51", 0xFF828282);}
		{int[] tmp = {133,133,133};
		colName2RGBdec.put("GREY52", tmp);
		colName2RGBhex.put("GREY52", 0xFF858585);}
		{int[] tmp = {135,135,135};
		colName2RGBdec.put("GREY53", tmp);
		colName2RGBhex.put("GREY53", 0xFF878787);}
		{int[] tmp = {138,138,138};
		colName2RGBdec.put("GREY54", tmp);
		colName2RGBhex.put("GREY54", 0xFF8A8A8A);}
		{int[] tmp = {140,140,140};
		colName2RGBdec.put("GREY55", tmp);
		colName2RGBhex.put("GREY55", 0xFF8C8C8C);}
		{int[] tmp = {143,143,143};
		colName2RGBdec.put("GREY56", tmp);
		colName2RGBhex.put("GREY56", 0xFF8F8F8F);}
		{int[] tmp = {145,145,145};
		colName2RGBdec.put("GREY57", tmp);
		colName2RGBhex.put("GREY57", 0xFF919191);}
		{int[] tmp = {148,148,148};
		colName2RGBdec.put("GREY58", tmp);
		colName2RGBhex.put("GREY58", 0xFF949494);}
		{int[] tmp = {150,150,150};
		colName2RGBdec.put("GREY59", tmp);
		colName2RGBhex.put("GREY59", 0xFF969696);}
		{int[] tmp = {153,153,153};
		colName2RGBdec.put("GREY60", tmp);
		colName2RGBhex.put("GREY60", 0xFF999999);}
		{int[] tmp = {156,156,156};
		colName2RGBdec.put("GREY61", tmp);
		colName2RGBhex.put("GREY61", 0xFF9C9C9C);}
		{int[] tmp = {158,158,158};
		colName2RGBdec.put("GREY62", tmp);
		colName2RGBhex.put("GREY62", 0xFF9E9E9E);}
		{int[] tmp = {161,161,161};
		colName2RGBdec.put("GREY63", tmp);
		colName2RGBhex.put("GREY63", 0xFFA1A1A1);}
		{int[] tmp = {163,163,163};
		colName2RGBdec.put("GREY64", tmp);
		colName2RGBhex.put("GREY64", 0xFFA3A3A3);}
		{int[] tmp = {166,166,166};
		colName2RGBdec.put("GREY65", tmp);
		colName2RGBhex.put("GREY65", 0xFFA6A6A6);}
		{int[] tmp = {168,168,168};
		colName2RGBdec.put("GREY66", tmp);
		colName2RGBhex.put("GREY66", 0xFFA8A8A8);}
		{int[] tmp = {171,171,171};
		colName2RGBdec.put("GREY67", tmp);
		colName2RGBhex.put("GREY67", 0xFFABABAB);}
		{int[] tmp = {173,173,173};
		colName2RGBdec.put("GREY68", tmp);
		colName2RGBhex.put("GREY68", 0xFFADADAD);}
		{int[] tmp = {176,176,176};
		colName2RGBdec.put("GREY69", tmp);
		colName2RGBhex.put("GREY69", 0xFFB0B0B0);}
		{int[] tmp = {179,179,179};
		colName2RGBdec.put("GREY70", tmp);
		colName2RGBhex.put("GREY70", 0xFFB3B3B3);}
		{int[] tmp = {181,181,181};
		colName2RGBdec.put("GREY71", tmp);
		colName2RGBhex.put("GREY71", 0xFFB5B5B5);}
		{int[] tmp = {184,184,184};
		colName2RGBdec.put("GREY72", tmp);
		colName2RGBhex.put("GREY72", 0xFFB8B8B8);}
		{int[] tmp = {186,186,186};
		colName2RGBdec.put("GREY73", tmp);
		colName2RGBhex.put("GREY73", 0xFFBABABA);}
		{int[] tmp = {189,189,189};
		colName2RGBdec.put("GREY74", tmp);
		colName2RGBhex.put("GREY74", 0xFFBDBDBD);}
		{int[] tmp = {191,191,191};
		colName2RGBdec.put("GREY75", tmp);
		colName2RGBhex.put("GREY75", 0xFFBFBFBF);}
		{int[] tmp = {194,194,194};
		colName2RGBdec.put("GREY76", tmp);
		colName2RGBhex.put("GREY76", 0xFFC2C2C2);}
		{int[] tmp = {196,196,196};
		colName2RGBdec.put("GREY77", tmp);
		colName2RGBhex.put("GREY77", 0xFFC4C4C4);}
		{int[] tmp = {199,199,199};
		colName2RGBdec.put("GREY78", tmp);
		colName2RGBhex.put("GREY78", 0xFFC7C7C7);}
		{int[] tmp = {201,201,201};
		colName2RGBdec.put("GREY79", tmp);
		colName2RGBhex.put("GREY79", 0xFFC9C9C9);}
		{int[] tmp = {204,204,204};
		colName2RGBdec.put("GREY80", tmp);
		colName2RGBhex.put("GREY80", 0xFFCCCCCC);}
		{int[] tmp = {207,207,207};
		colName2RGBdec.put("GREY81", tmp);
		colName2RGBhex.put("GREY81", 0xFFCFCFCF);}
		{int[] tmp = {209,209,209};
		colName2RGBdec.put("GREY82", tmp);
		colName2RGBhex.put("GREY82", 0xFFD1D1D1);}
		{int[] tmp = {212,212,212};
		colName2RGBdec.put("GREY83", tmp);
		colName2RGBhex.put("GREY83", 0xFFD4D4D4);}
		{int[] tmp = {214,214,214};
		colName2RGBdec.put("GREY84", tmp);
		colName2RGBhex.put("GREY84", 0xFFD6D6D6);}
		{int[] tmp = {217,217,217};
		colName2RGBdec.put("GREY85", tmp);
		colName2RGBhex.put("GREY85", 0xFFD9D9D9);}
		{int[] tmp = {219,219,219};
		colName2RGBdec.put("GREY86", tmp);
		colName2RGBhex.put("GREY86", 0xFFDBDBDB);}
		{int[] tmp = {222,222,222};
		colName2RGBdec.put("GREY87", tmp);
		colName2RGBhex.put("GREY87", 0xFFDEDEDE);}
		{int[] tmp = {224,224,224};
		colName2RGBdec.put("GREY88", tmp);
		colName2RGBhex.put("GREY88", 0xFFE0E0E0);}
		{int[] tmp = {227,227,227};
		colName2RGBdec.put("GREY89", tmp);
		colName2RGBhex.put("GREY89", 0xFFE3E3E3);}
		{int[] tmp = {229,229,229};
		colName2RGBdec.put("GREY90", tmp);
		colName2RGBhex.put("GREY90", 0xFFE5E5E5);}
		{int[] tmp = {232,232,232};
		colName2RGBdec.put("GREY91", tmp);
		colName2RGBhex.put("GREY91", 0xFFE8E8E8);}
		{int[] tmp = {235,235,235};
		colName2RGBdec.put("GREY92", tmp);
		colName2RGBhex.put("GREY92", 0xFFEBEBEB);}
		{int[] tmp = {237,237,237};
		colName2RGBdec.put("GREY93", tmp);
		colName2RGBhex.put("GREY93", 0xFFEDEDED);}
		{int[] tmp = {240,240,240};
		colName2RGBdec.put("GREY94", tmp);
		colName2RGBhex.put("GREY94", 0xFFF0F0F0);}
		{int[] tmp = {242,242,242};
		colName2RGBdec.put("GREY95", tmp);
		colName2RGBhex.put("GREY95", 0xFFF2F2F2);}
		{int[] tmp = {245,245,245};
		colName2RGBdec.put("GREY96", tmp);
		colName2RGBhex.put("GREY96", 0xFFF5F5F5);}
		{int[] tmp = {247,247,247};
		colName2RGBdec.put("GREY97", tmp);
		colName2RGBhex.put("GREY97", 0xFFF7F7F7);}
		{int[] tmp = {250,250,250};
		colName2RGBdec.put("GREY98", tmp);
		colName2RGBhex.put("GREY98", 0xFFFAFAFA);}
		{int[] tmp = {252,252,252};
		colName2RGBdec.put("GREY99", tmp);
		colName2RGBhex.put("GREY99", 0xFFFCFCFC);}
		{int[] tmp = {255,255,255};
		colName2RGBdec.put("GREY100", tmp);
		colName2RGBhex.put("GREY100", 0xFFFFFFFF);}
		{int[] tmp = {255,255,255};
		colName2RGBdec.put("WHITE", tmp);
		colName2RGBhex.put("WHITE", 0xFFFFFFFF);}
		{int[] tmp = {240,248,255};
		colName2RGBdec.put("ALICEBLUE", tmp);
		colName2RGBhex.put("ALICEBLUE", 0xFFF0F8FF);}
		{int[] tmp = {138,43,226};
		colName2RGBdec.put("BLUEVIOLET", tmp);
		colName2RGBhex.put("BLUEVIOLET", 0xFF8A2BE2);}
		{int[] tmp = {95,158,160};
		colName2RGBdec.put("CADETBLUE", tmp);
		colName2RGBhex.put("CADETBLUE", 0xFF5F9EA0);}
		{int[] tmp = {152,245,255};
		colName2RGBdec.put("CADETBLUE1", tmp);
		colName2RGBhex.put("CADETBLUE1", 0xFF98F5FF);}
		{int[] tmp = {142,229,238};
		colName2RGBdec.put("CADETBLUE2", tmp);
		colName2RGBhex.put("CADETBLUE2", 0xFF8EE5EE);}
		{int[] tmp = {122,197,205};
		colName2RGBdec.put("CADETBLUE3", tmp);
		colName2RGBhex.put("CADETBLUE3", 0xFF7AC5CD);}
		{int[] tmp = {83,134,139};
		colName2RGBdec.put("CADETBLUE4", tmp);
		colName2RGBhex.put("CADETBLUE4", 0xFF53868B);}
		{int[] tmp = {100,149,237};
		colName2RGBdec.put("CORNFLOWERBLUE", tmp);
		colName2RGBhex.put("CORNFLOWERBLUE", 0xFF6495ED);}
		{int[] tmp = {72,61,139};
		colName2RGBdec.put("DARKSLATEBLUE", tmp);
		colName2RGBhex.put("DARKSLATEBLUE", 0xFF483D8B);}
		{int[] tmp = {0,206,209};
		colName2RGBdec.put("DARKTURQUOISE", tmp);
		colName2RGBhex.put("DARKTURQUOISE", 0xFF00CED1);}
		{int[] tmp = {0,191,255};
		colName2RGBdec.put("DEEPSKYBLUE", tmp);
		colName2RGBhex.put("DEEPSKYBLUE", 0xFF00BFFF);}
		{int[] tmp = {0,191,255};
		colName2RGBdec.put("DEEPSKYBLUE1", tmp);
		colName2RGBhex.put("DEEPSKYBLUE1", 0xFF00BFFF);}
		{int[] tmp = {0,178,238};
		colName2RGBdec.put("DEEPSKYBLUE2", tmp);
		colName2RGBhex.put("DEEPSKYBLUE2", 0xFF00B2EE);}
		{int[] tmp = {0,154,205};
		colName2RGBdec.put("DEEPSKYBLUE3", tmp);
		colName2RGBhex.put("DEEPSKYBLUE3", 0xFF009ACD);}
		{int[] tmp = {0,104,139};
		colName2RGBdec.put("DEEPSKYBLUE4", tmp);
		colName2RGBhex.put("DEEPSKYBLUE4", 0xFF00688B);}
		{int[] tmp = {30,144,255};
		colName2RGBdec.put("DODGERBLUE", tmp);
		colName2RGBhex.put("DODGERBLUE", 0xFF1E90FF);}
		{int[] tmp = {30,144,255};
		colName2RGBdec.put("DODGERBLUE1", tmp);
		colName2RGBhex.put("DODGERBLUE1", 0xFF1E90FF);}
		{int[] tmp = {28,134,238};
		colName2RGBdec.put("DODGERBLUE2", tmp);
		colName2RGBhex.put("DODGERBLUE2", 0xFF1C86EE);}
		{int[] tmp = {24,116,205};
		colName2RGBdec.put("DODGERBLUE3", tmp);
		colName2RGBhex.put("DODGERBLUE3", 0xFF1874CD);}
		{int[] tmp = {16,78,139};
		colName2RGBdec.put("DODGERBLUE4", tmp);
		colName2RGBhex.put("DODGERBLUE4", 0xFF104E8B);}
		{int[] tmp = {173,216,230};
		colName2RGBdec.put("LIGHTBLUE", tmp);
		colName2RGBhex.put("LIGHTBLUE", 0xFFADD8E6);}
		{int[] tmp = {191,239,255};
		colName2RGBdec.put("LIGHTBLUE1", tmp);
		colName2RGBhex.put("LIGHTBLUE1", 0xFFBFEFFF);}
		{int[] tmp = {178,223,238};
		colName2RGBdec.put("LIGHTBLUE2", tmp);
		colName2RGBhex.put("LIGHTBLUE2", 0xFFB2DFEE);}
		{int[] tmp = {154,192,205};
		colName2RGBdec.put("LIGHTBLUE3", tmp);
		colName2RGBhex.put("LIGHTBLUE3", 0xFF9AC0CD);}
		{int[] tmp = {104,131,139};
		colName2RGBdec.put("LIGHTBLUE4", tmp);
		colName2RGBhex.put("LIGHTBLUE4", 0xFF68838B);}
		{int[] tmp = {224,255,255};
		colName2RGBdec.put("LIGHTCYAN", tmp);
		colName2RGBhex.put("LIGHTCYAN", 0xFFE0FFFF);}
		{int[] tmp = {224,255,255};
		colName2RGBdec.put("LIGHTCYAN1", tmp);
		colName2RGBhex.put("LIGHTCYAN1", 0xFFE0FFFF);}
		{int[] tmp = {209,238,238};
		colName2RGBdec.put("LIGHTCYAN2", tmp);
		colName2RGBhex.put("LIGHTCYAN2", 0xFFD1EEEE);}
		{int[] tmp = {180,205,205};
		colName2RGBdec.put("LIGHTCYAN3", tmp);
		colName2RGBhex.put("LIGHTCYAN3", 0xFFB4CDCD);}
		{int[] tmp = {122,139,139};
		colName2RGBdec.put("LIGHTCYAN4", tmp);
		colName2RGBhex.put("LIGHTCYAN4", 0xFF7A8B8B);}
		{int[] tmp = {135,206,250};
		colName2RGBdec.put("LIGHTSKYBLUE", tmp);
		colName2RGBhex.put("LIGHTSKYBLUE", 0xFF87CEFA);}
		{int[] tmp = {176,226,255};
		colName2RGBdec.put("LIGHTSKYBLUE1", tmp);
		colName2RGBhex.put("LIGHTSKYBLUE1", 0xFFB0E2FF);}
		{int[] tmp = {164,211,238};
		colName2RGBdec.put("LIGHTSKYBLUE2", tmp);
		colName2RGBhex.put("LIGHTSKYBLUE2", 0xFFA4D3EE);}
		{int[] tmp = {141,182,205};
		colName2RGBdec.put("LIGHTSKYBLUE3", tmp);
		colName2RGBhex.put("LIGHTSKYBLUE3", 0xFF8DB6CD);}
		{int[] tmp = {96,123,139};
		colName2RGBdec.put("LIGHTSKYBLUE4", tmp);
		colName2RGBhex.put("LIGHTSKYBLUE4", 0xFF607B8B);}
		{int[] tmp = {132,112,255};
		colName2RGBdec.put("LIGHTSLATEBLUE", tmp);
		colName2RGBhex.put("LIGHTSLATEBLUE", 0xFF8470FF);}
		{int[] tmp = {176,196,222};
		colName2RGBdec.put("LIGHTSTEELBLUE", tmp);
		colName2RGBhex.put("LIGHTSTEELBLUE", 0xFFB0C4DE);}
		{int[] tmp = {202,225,255};
		colName2RGBdec.put("LIGHTSTEELBLUE1", tmp);
		colName2RGBhex.put("LIGHTSTEELBLUE1", 0xFFCAE1FF);}
		{int[] tmp = {188,210,238};
		colName2RGBdec.put("LIGHTSTEELBLUE2", tmp);
		colName2RGBhex.put("LIGHTSTEELBLUE2", 0xFFBCD2EE);}
		{int[] tmp = {162,181,205};
		colName2RGBdec.put("LIGHTSTEELBLUE3", tmp);
		colName2RGBhex.put("LIGHTSTEELBLUE3", 0xFFA2B5CD);}
		{int[] tmp = {110,123,139};
		colName2RGBdec.put("LIGHTSTEELBLUE4", tmp);
		colName2RGBhex.put("LIGHTSTEELBLUE4", 0xFF6E7B8B);}
		{int[] tmp = {0,0,205};
		colName2RGBdec.put("MEDIUMBLUE", tmp);
		colName2RGBhex.put("MEDIUMBLUE", 0xFF0000CD);}
		{int[] tmp = {123,104,238};
		colName2RGBdec.put("MEDIUMSLATEBLUE", tmp);
		colName2RGBhex.put("MEDIUMSLATEBLUE", 0xFF7B68EE);}
		{int[] tmp = {72,209,204};
		colName2RGBdec.put("MEDIUMTURQUOISE", tmp);
		colName2RGBhex.put("MEDIUMTURQUOISE", 0xFF48D1CC);}
		{int[] tmp = {25,25,112};
		colName2RGBdec.put("MIDNIGHTBLUE", tmp);
		colName2RGBhex.put("MIDNIGHTBLUE", 0xFF191970);}
		{int[] tmp = {0,0,128};
		colName2RGBdec.put("NAVYBLUE", tmp);
		colName2RGBhex.put("NAVYBLUE", 0xFF000080);}
		{int[] tmp = {175,238,238};
		colName2RGBdec.put("PALETURQUOISE", tmp);
		colName2RGBhex.put("PALETURQUOISE", 0xFFAFEEEE);}
		{int[] tmp = {187,255,255};
		colName2RGBdec.put("PALETURQUOISE1", tmp);
		colName2RGBhex.put("PALETURQUOISE1", 0xFFBBFFFF);}
		{int[] tmp = {174,238,238};
		colName2RGBdec.put("PALETURQUOISE2", tmp);
		colName2RGBhex.put("PALETURQUOISE2", 0xFFAEEEEE);}
		{int[] tmp = {150,205,205};
		colName2RGBdec.put("PALETURQUOISE3", tmp);
		colName2RGBhex.put("PALETURQUOISE3", 0xFF96CDCD);}
		{int[] tmp = {102,139,139};
		colName2RGBdec.put("PALETURQUOISE4", tmp);
		colName2RGBhex.put("PALETURQUOISE4", 0xFF668B8B);}
		{int[] tmp = {176,224,230};
		colName2RGBdec.put("POWDERBLUE", tmp);
		colName2RGBhex.put("POWDERBLUE", 0xFFB0E0E6);}
		{int[] tmp = {65,105,225};
		colName2RGBdec.put("ROYALBLUE", tmp);
		colName2RGBhex.put("ROYALBLUE", 0xFF4169E1);}
		{int[] tmp = {72,118,255};
		colName2RGBdec.put("ROYALBLUE1", tmp);
		colName2RGBhex.put("ROYALBLUE1", 0xFF4876FF);}
		{int[] tmp = {67,110,238};
		colName2RGBdec.put("ROYALBLUE2", tmp);
		colName2RGBhex.put("ROYALBLUE2", 0xFF436EEE);}
		{int[] tmp = {58,95,205};
		colName2RGBdec.put("ROYALBLUE3", tmp);
		colName2RGBhex.put("ROYALBLUE3", 0xFF3A5FCD);}
		{int[] tmp = {39,64,139};
		colName2RGBdec.put("ROYALBLUE4", tmp);
		colName2RGBhex.put("ROYALBLUE4", 0xFF27408B);}
		{int[] tmp = {0,34,102};
		colName2RGBdec.put("ROYALBLUE5", tmp);
		colName2RGBhex.put("ROYALBLUE5", 0xFF002266);}
		{int[] tmp = {135,206,235};
		colName2RGBdec.put("SKYBLUE", tmp);
		colName2RGBhex.put("SKYBLUE", 0xFF87CEEB);}
		{int[] tmp = {135,206,255};
		colName2RGBdec.put("SKYBLUE1", tmp);
		colName2RGBhex.put("SKYBLUE1", 0xFF87CEFF);}
		{int[] tmp = {126,192,238};
		colName2RGBdec.put("SKYBLUE2", tmp);
		colName2RGBhex.put("SKYBLUE2", 0xFF7EC0EE);}
		{int[] tmp = {108,166,205};
		colName2RGBdec.put("SKYBLUE3", tmp);
		colName2RGBhex.put("SKYBLUE3", 0xFF6CA6CD);}
		{int[] tmp = {74,112,139};
		colName2RGBdec.put("SKYBLUE4", tmp);
		colName2RGBhex.put("SKYBLUE4", 0xFF4A708B);}
		{int[] tmp = {106,90,205};
		colName2RGBdec.put("SLATEBLUE", tmp);
		colName2RGBhex.put("SLATEBLUE", 0xFF6A5ACD);}
		{int[] tmp = {131,111,255};
		colName2RGBdec.put("SLATEBLUE1", tmp);
		colName2RGBhex.put("SLATEBLUE1", 0xFF836FFF);}
		{int[] tmp = {122,103,238};
		colName2RGBdec.put("SLATEBLUE2", tmp);
		colName2RGBhex.put("SLATEBLUE2", 0xFF7A67EE);}
		{int[] tmp = {105,89,205};
		colName2RGBdec.put("SLATEBLUE3", tmp);
		colName2RGBhex.put("SLATEBLUE3", 0xFF6959CD);}
		{int[] tmp = {71,60,139};
		colName2RGBdec.put("SLATEBLUE4", tmp);
		colName2RGBhex.put("SLATEBLUE4", 0xFF473C8B);}
		{int[] tmp = {70,130,180};
		colName2RGBdec.put("STEELBLUE", tmp);
		colName2RGBhex.put("STEELBLUE", 0xFF4682B4);}
		{int[] tmp = {99,184,255};
		colName2RGBdec.put("STEELBLUE1", tmp);
		colName2RGBhex.put("STEELBLUE1", 0xFF63B8FF);}
		{int[] tmp = {92,172,238};
		colName2RGBdec.put("STEELBLUE2", tmp);
		colName2RGBhex.put("STEELBLUE2", 0xFF5CACEE);}
		{int[] tmp = {79,148,205};
		colName2RGBdec.put("STEELBLUE3", tmp);
		colName2RGBhex.put("STEELBLUE3", 0xFF4F94CD);}
		{int[] tmp = {54,100,139};
		colName2RGBdec.put("STEELBLUE4", tmp);
		colName2RGBhex.put("STEELBLUE4", 0xFF36648B);}
		{int[] tmp = {127,255,212};
		colName2RGBdec.put("AQUAMARINE", tmp);
		colName2RGBhex.put("AQUAMARINE", 0xFF7FFFD4);}
		{int[] tmp = {127,255,212};
		colName2RGBdec.put("AQUAMARINE1", tmp);
		colName2RGBhex.put("AQUAMARINE1", 0xFF7FFFD4);}
		{int[] tmp = {118,238,198};
		colName2RGBdec.put("AQUAMARINE2", tmp);
		colName2RGBhex.put("AQUAMARINE2", 0xFF76EEC6);}
		{int[] tmp = {102,205,170};
		colName2RGBdec.put("AQUAMARINE3", tmp);
		colName2RGBhex.put("AQUAMARINE3", 0xFF66CDAA);}
		{int[] tmp = {102,205,170};
		colName2RGBdec.put("MEDIUMAQUAMARINE", tmp);
		colName2RGBhex.put("MEDIUMAQUAMARINE", 0xFF66CDAA);}
		{int[] tmp = {69,139,116};
		colName2RGBdec.put("AQUAMARINE4", tmp);
		colName2RGBhex.put("AQUAMARINE4", 0xFF458B74);}
		{int[] tmp = {240,255,255};
		colName2RGBdec.put("AZURE", tmp);
		colName2RGBhex.put("AZURE", 0xFFF0FFFF);}
		{int[] tmp = {240,255,255};
		colName2RGBdec.put("AZURE1", tmp);
		colName2RGBhex.put("AZURE1", 0xFFF0FFFF);}
		{int[] tmp = {224,238,238};
		colName2RGBdec.put("AZURE2", tmp);
		colName2RGBhex.put("AZURE2", 0xFFE0EEEE);}
		{int[] tmp = {193,205,205};
		colName2RGBdec.put("AZURE3", tmp);
		colName2RGBhex.put("AZURE3", 0xFFC1CDCD);}
		{int[] tmp = {131,139,139};
		colName2RGBdec.put("AZURE4", tmp);
		colName2RGBhex.put("AZURE4", 0xFF838B8B);}
		{int[] tmp = {0,0,255};
		colName2RGBdec.put("BLUE", tmp);
		colName2RGBhex.put("BLUE", 0xFF0000FF);}
		{int[] tmp = {0,0,255};
		colName2RGBdec.put("BLUE1", tmp);
		colName2RGBhex.put("BLUE1", 0xFF0000FF);}
		{int[] tmp = {0,0,238};
		colName2RGBdec.put("BLUE2", tmp);
		colName2RGBhex.put("BLUE2", 0xFF0000EE);}
		{int[] tmp = {0,0,205};
		colName2RGBdec.put("BLUE3", tmp);
		colName2RGBhex.put("BLUE3", 0xFF0000CD);}
		{int[] tmp = {0,0,139};
		colName2RGBdec.put("BLUE4", tmp);
		colName2RGBhex.put("BLUE4", 0xFF00008B);}
		{int[] tmp = {0,255,255};
		colName2RGBdec.put("CYAN", tmp);
		colName2RGBhex.put("CYAN", 0xFF00FFFF);}
		{int[] tmp = {0,255,255};
		colName2RGBdec.put("CYAN1", tmp);
		colName2RGBhex.put("CYAN1", 0xFF00FFFF);}
		{int[] tmp = {0,238,238};
		colName2RGBdec.put("CYAN2", tmp);
		colName2RGBhex.put("CYAN2", 0xFF00EEEE);}
		{int[] tmp = {0,205,205};
		colName2RGBdec.put("CYAN3", tmp);
		colName2RGBhex.put("CYAN3", 0xFF00CDCD);}
		{int[] tmp = {0,139,139};
		colName2RGBdec.put("CYAN4", tmp);
		colName2RGBhex.put("CYAN4", 0xFF008B8B);}
		{int[] tmp = {0,0,128};
		colName2RGBdec.put("NAVY", tmp);
		colName2RGBhex.put("NAVY", 0xFF000080);}
		{int[] tmp = {64,224,208};
		colName2RGBdec.put("TURQUOISE", tmp);
		colName2RGBhex.put("TURQUOISE", 0xFF40E0D0);}
		{int[] tmp = {0,245,255};
		colName2RGBdec.put("TURQUOISE1", tmp);
		colName2RGBhex.put("TURQUOISE1", 0xFF00F5FF);}
		{int[] tmp = {0,229,238};
		colName2RGBdec.put("TURQUOISE2", tmp);
		colName2RGBhex.put("TURQUOISE2", 0xFF00E5EE);}
		{int[] tmp = {0,197,205};
		colName2RGBdec.put("TURQUOISE3", tmp);
		colName2RGBhex.put("TURQUOISE3", 0xFF00C5CD);}
		{int[] tmp = {0,134,139};
		colName2RGBdec.put("TURQUOISE4", tmp);
		colName2RGBhex.put("TURQUOISE4", 0xFF00868B);}
		{int[] tmp = {47,79,79};
		colName2RGBdec.put("DARKSLATEGRAY", tmp);
		colName2RGBhex.put("DARKSLATEGRAY", 0xFF2F4F4F);}
		{int[] tmp = {151,255,255};
		colName2RGBdec.put("DARKSLATEGRAY1", tmp);
		colName2RGBhex.put("DARKSLATEGRAY1", 0xFF97FFFF);}
		{int[] tmp = {141,238,238};
		colName2RGBdec.put("DARKSLATEGRAY2", tmp);
		colName2RGBhex.put("DARKSLATEGRAY2", 0xFF8DEEEE);}
		{int[] tmp = {121,205,205};
		colName2RGBdec.put("DARKSLATEGRAY3", tmp);
		colName2RGBhex.put("DARKSLATEGRAY3", 0xFF79CDCD);}
		{int[] tmp = {82,139,139};
		colName2RGBdec.put("DARKSLATEGRAY4", tmp);
		colName2RGBhex.put("DARKSLATEGRAY4", 0xFF528B8B);}
		{int[] tmp = {188,143,143};
		colName2RGBdec.put("ROSYBROWN", tmp);
		colName2RGBhex.put("ROSYBROWN", 0xFFBC8F8F);}
		{int[] tmp = {255,193,193};
		colName2RGBdec.put("ROSYBROWN1", tmp);
		colName2RGBhex.put("ROSYBROWN1", 0xFFFFC1C1);}
		{int[] tmp = {238,180,180};
		colName2RGBdec.put("ROSYBROWN2", tmp);
		colName2RGBhex.put("ROSYBROWN2", 0xFFEEB4B4);}
		{int[] tmp = {205,155,155};
		colName2RGBdec.put("ROSYBROWN3", tmp);
		colName2RGBhex.put("ROSYBROWN3", 0xFFCD9B9B);}
		{int[] tmp = {139,105,105};
		colName2RGBdec.put("ROSYBROWN4", tmp);
		colName2RGBhex.put("ROSYBROWN4", 0xFF8B6969);}
		{int[] tmp = {139,69,19};
		colName2RGBdec.put("SADDLEBROWN", tmp);
		colName2RGBhex.put("SADDLEBROWN", 0xFF8B4513);}
		{int[] tmp = {244,164,96};
		colName2RGBdec.put("SANDYBROWN", tmp);
		colName2RGBhex.put("SANDYBROWN", 0xFFF4A460);}
		{int[] tmp = {245,245,220};
		colName2RGBdec.put("BEIGE", tmp);
		colName2RGBhex.put("BEIGE", 0xFFF5F5DC);}
		{int[] tmp = {165,42,42};
		colName2RGBdec.put("BROWN", tmp);
		colName2RGBhex.put("BROWN", 0xFFA52A2A);}
		{int[] tmp = {255,64,64};
		colName2RGBdec.put("BROWN1", tmp);
		colName2RGBhex.put("BROWN1", 0xFFFF4040);}
		{int[] tmp = {238,59,59};
		colName2RGBdec.put("BROWN2", tmp);
		colName2RGBhex.put("BROWN2", 0xFFEE3B3B);}
		{int[] tmp = {205,51,51};
		colName2RGBdec.put("BROWN3", tmp);
		colName2RGBhex.put("BROWN3", 0xFFCD3333);}
		{int[] tmp = {139,35,35};
		colName2RGBdec.put("BROWN4", tmp);
		colName2RGBhex.put("BROWN4", 0xFF8B2323);}
		{int[] tmp = {222,184,135};
		colName2RGBdec.put("BURLYWOOD", tmp);
		colName2RGBhex.put("BURLYWOOD", 0xFFDEB887);}
		{int[] tmp = {255,211,155};
		colName2RGBdec.put("BURLYWOOD1", tmp);
		colName2RGBhex.put("BURLYWOOD1", 0xFFFFD39B);}
		{int[] tmp = {238,197,145};
		colName2RGBdec.put("BURLYWOOD2", tmp);
		colName2RGBhex.put("BURLYWOOD2", 0xFFEEC591);}
		{int[] tmp = {205,170,125};
		colName2RGBdec.put("BURLYWOOD3", tmp);
		colName2RGBhex.put("BURLYWOOD3", 0xFFCDAA7D);}
		{int[] tmp = {139,115,85};
		colName2RGBdec.put("BURLYWOOD4", tmp);
		colName2RGBhex.put("BURLYWOOD4", 0xFF8B7355);}
		{int[] tmp = {210,105,30};
		colName2RGBdec.put("CHOCOLATE", tmp);
		colName2RGBhex.put("CHOCOLATE", 0xFFD2691E);}
		{int[] tmp = {255,127,36};
		colName2RGBdec.put("CHOCOLATE1", tmp);
		colName2RGBhex.put("CHOCOLATE1", 0xFFFF7F24);}
		{int[] tmp = {238,118,33};
		colName2RGBdec.put("CHOCOLATE2", tmp);
		colName2RGBhex.put("CHOCOLATE2", 0xFFEE7621);}
		{int[] tmp = {205,102,29};
		colName2RGBdec.put("CHOCOLATE3", tmp);
		colName2RGBhex.put("CHOCOLATE3", 0xFFCD661D);}
		{int[] tmp = {139,69,19};
		colName2RGBdec.put("CHOCOLATE4", tmp);
		colName2RGBhex.put("CHOCOLATE4", 0xFF8B4513);}
		{int[] tmp = {205,133,63};
		colName2RGBdec.put("PERU", tmp);
		colName2RGBhex.put("PERU", 0xFFCD853F);}
		{int[] tmp = {210,180,140};
		colName2RGBdec.put("TAN", tmp);
		colName2RGBhex.put("TAN", 0xFFD2B48C);}
		{int[] tmp = {255,165,79};
		colName2RGBdec.put("TAN1", tmp);
		colName2RGBhex.put("TAN1", 0xFFFFA54F);}
		{int[] tmp = {238,154,73};
		colName2RGBdec.put("TAN2", tmp);
		colName2RGBhex.put("TAN2", 0xFFEE9A49);}
		{int[] tmp = {205,133,63};
		colName2RGBdec.put("TAN3", tmp);
		colName2RGBhex.put("TAN3", 0xFFCD853F);}
		{int[] tmp = {139,90,43};
		colName2RGBdec.put("TAN4", tmp);
		colName2RGBhex.put("TAN4", 0xFF8B5A2B);}
		{int[] tmp = {0,100,0};
		colName2RGBdec.put("DARKGREEN", tmp);
		colName2RGBhex.put("DARKGREEN", 0xFF006400);}
		{int[] tmp = {189,183,107};
		colName2RGBdec.put("DARKKHAKI", tmp);
		colName2RGBhex.put("DARKKHAKI", 0xFFBDB76B);}
		{int[] tmp = {85,107,47};
		colName2RGBdec.put("DARKOLIVEGREEN", tmp);
		colName2RGBhex.put("DARKOLIVEGREEN", 0xFF556B2F);}
		{int[] tmp = {202,255,112};
		colName2RGBdec.put("DARKOLIVEGREEN1", tmp);
		colName2RGBhex.put("DARKOLIVEGREEN1", 0xFFCAFF70);}
		{int[] tmp = {188,238,104};
		colName2RGBdec.put("DARKOLIVEGREEN2", tmp);
		colName2RGBhex.put("DARKOLIVEGREEN2", 0xFFBCEE68);}
		{int[] tmp = {162,205,90};
		colName2RGBdec.put("DARKOLIVEGREEN3", tmp);
		colName2RGBhex.put("DARKOLIVEGREEN3", 0xFFA2CD5A);}
		{int[] tmp = {110,139,61};
		colName2RGBdec.put("DARKOLIVEGREEN4", tmp);
		colName2RGBhex.put("DARKOLIVEGREEN4", 0xFF6E8B3D);}
		{int[] tmp = {143,188,143};
		colName2RGBdec.put("DARKSEAGREEN", tmp);
		colName2RGBhex.put("DARKSEAGREEN", 0xFF8FBC8F);}
		{int[] tmp = {193,255,193};
		colName2RGBdec.put("DARKSEAGREEN1", tmp);
		colName2RGBhex.put("DARKSEAGREEN1", 0xFFC1FFC1);}
		{int[] tmp = {180,238,180};
		colName2RGBdec.put("DARKSEAGREEN2", tmp);
		colName2RGBhex.put("DARKSEAGREEN2", 0xFFB4EEB4);}
		{int[] tmp = {155,205,155};
		colName2RGBdec.put("DARKSEAGREEN3", tmp);
		colName2RGBhex.put("DARKSEAGREEN3", 0xFF9BCD9B);}
		{int[] tmp = {105,139,105};
		colName2RGBdec.put("DARKSEAGREEN4", tmp);
		colName2RGBhex.put("DARKSEAGREEN4", 0xFF698B69);}
		{int[] tmp = {34,139,34};
		colName2RGBdec.put("FORESTGREEN", tmp);
		colName2RGBhex.put("FORESTGREEN", 0xFF228B22);}
		{int[] tmp = {173,255,47};
		colName2RGBdec.put("GREENYELLOW", tmp);
		colName2RGBhex.put("GREENYELLOW", 0xFFADFF2F);}
		{int[] tmp = {124,252,0};
		colName2RGBdec.put("LAWNGREEN", tmp);
		colName2RGBhex.put("LAWNGREEN", 0xFF7CFC00);}
		{int[] tmp = {32,178,170};
		colName2RGBdec.put("LIGHTSEAGREEN", tmp);
		colName2RGBhex.put("LIGHTSEAGREEN", 0xFF20B2AA);}
		{int[] tmp = {50,205,50};
		colName2RGBdec.put("LIMEGREEN", tmp);
		colName2RGBhex.put("LIMEGREEN", 0xFF32CD32);}
		{int[] tmp = {60,179,113};
		colName2RGBdec.put("MEDIUMSEAGREEN", tmp);
		colName2RGBhex.put("MEDIUMSEAGREEN", 0xFF3CB371);}
		{int[] tmp = {0,250,154};
		colName2RGBdec.put("MEDIUMSPRINGGREEN", tmp);
		colName2RGBhex.put("MEDIUMSPRINGGREEN", 0xFF00FA9A);}
		{int[] tmp = {245,255,250};
		colName2RGBdec.put("MINTCREAM", tmp);
		colName2RGBhex.put("MINTCREAM", 0xFFF5FFFA);}
		{int[] tmp = {107,142,35};
		colName2RGBdec.put("OLIVEDRAB", tmp);
		colName2RGBhex.put("OLIVEDRAB", 0xFF6B8E23);}
		{int[] tmp = {192,255,62};
		colName2RGBdec.put("OLIVEDRAB1", tmp);
		colName2RGBhex.put("OLIVEDRAB1", 0xFFC0FF3E);}
		{int[] tmp = {179,238,58};
		colName2RGBdec.put("OLIVEDRAB2", tmp);
		colName2RGBhex.put("OLIVEDRAB2", 0xFFB3EE3A);}
		{int[] tmp = {154,205,50};
		colName2RGBdec.put("OLIVEDRAB3", tmp);
		colName2RGBhex.put("OLIVEDRAB3", 0xFF9ACD32);}
		{int[] tmp = {105,139,34};
		colName2RGBdec.put("OLIVEDRAB4", tmp);
		colName2RGBhex.put("OLIVEDRAB4", 0xFF698B22);}
		{int[] tmp = {152,251,152};
		colName2RGBdec.put("PALEGREEN", tmp);
		colName2RGBhex.put("PALEGREEN", 0xFF98FB98);}
		{int[] tmp = {154,255,154};
		colName2RGBdec.put("PALEGREEN1", tmp);
		colName2RGBhex.put("PALEGREEN1", 0xFF9AFF9A);}
		{int[] tmp = {144,238,144};
		colName2RGBdec.put("PALEGREEN2", tmp);
		colName2RGBhex.put("PALEGREEN2", 0xFF90EE90);}
		{int[] tmp = {124,205,124};
		colName2RGBdec.put("PALEGREEN3", tmp);
		colName2RGBhex.put("PALEGREEN3", 0xFF7CCD7C);}
		{int[] tmp = {84,139,84};
		colName2RGBdec.put("PALEGREEN4", tmp);
		colName2RGBhex.put("PALEGREEN4", 0xFF548B54);}
		{int[] tmp = {46,139,87};
		colName2RGBdec.put("SEAGREEN", tmp);
		colName2RGBhex.put("SEAGREEN", 0xFF2E8B57);}
		{int[] tmp = {46,139,87};
		colName2RGBdec.put("SEAGREEN4", tmp);
		colName2RGBhex.put("SEAGREEN4", 0xFF2E8B57);}
		{int[] tmp = {84,255,159};
		colName2RGBdec.put("SEAGREEN1", tmp);
		colName2RGBhex.put("SEAGREEN1", 0xFF54FF9F);}
		{int[] tmp = {78,238,148};
		colName2RGBdec.put("SEAGREEN2", tmp);
		colName2RGBhex.put("SEAGREEN2", 0xFF4EEE94);}
		{int[] tmp = {67,205,128};
		colName2RGBdec.put("SEAGREEN3", tmp);
		colName2RGBhex.put("SEAGREEN3", 0xFF43CD80);}
		{int[] tmp = {0,255,127};
		colName2RGBdec.put("SPRINGGREEN", tmp);
		colName2RGBhex.put("SPRINGGREEN", 0xFF00FF7F);}
		{int[] tmp = {0,255,127};
		colName2RGBdec.put("SPRINGGREEN1", tmp);
		colName2RGBhex.put("SPRINGGREEN1", 0xFF00FF7F);}
		{int[] tmp = {0,238,118};
		colName2RGBdec.put("SPRINGGREEN2", tmp);
		colName2RGBhex.put("SPRINGGREEN2", 0xFF00EE76);}
		{int[] tmp = {0,205,102};
		colName2RGBdec.put("SPRINGGREEN3", tmp);
		colName2RGBhex.put("SPRINGGREEN3", 0xFF00CD66);}
		{int[] tmp = {0,139,69};
		colName2RGBdec.put("SPRINGGREEN4", tmp);
		colName2RGBhex.put("SPRINGGREEN4", 0xFF008B45);}
		{int[] tmp = {154,205,50};
		colName2RGBdec.put("YELLOWGREEN", tmp);
		colName2RGBhex.put("YELLOWGREEN", 0xFF9ACD32);}
		{int[] tmp = {127,255,0};
		colName2RGBdec.put("CHARTREUSE", tmp);
		colName2RGBhex.put("CHARTREUSE", 0xFF7FFF00);}
		{int[] tmp = {127,255,0};
		colName2RGBdec.put("CHARTREUSE1", tmp);
		colName2RGBhex.put("CHARTREUSE1", 0xFF7FFF00);}
		{int[] tmp = {118,238,0};
		colName2RGBdec.put("CHARTREUSE2", tmp);
		colName2RGBhex.put("CHARTREUSE2", 0xFF76EE00);}
		{int[] tmp = {102,205,0};
		colName2RGBdec.put("CHARTREUSE3", tmp);
		colName2RGBhex.put("CHARTREUSE3", 0xFF66CD00);}
		{int[] tmp = {69,139,0};
		colName2RGBdec.put("CHARTREUSE4", tmp);
		colName2RGBhex.put("CHARTREUSE4", 0xFF458B00);}
		{int[] tmp = {0,255,0};
		colName2RGBdec.put("GREEN", tmp);
		colName2RGBhex.put("GREEN", 0xFF00FF00);}
		{int[] tmp = {0,255,0};
		colName2RGBdec.put("GREEN1", tmp);
		colName2RGBhex.put("GREEN1", 0xFF00FF00);}
		{int[] tmp = {0,238,0};
		colName2RGBdec.put("GREEN2", tmp);
		colName2RGBhex.put("GREEN2", 0xFF00EE00);}
		{int[] tmp = {0,205,0};
		colName2RGBdec.put("GREEN3", tmp);
		colName2RGBhex.put("GREEN3", 0xFF00CD00);}
		{int[] tmp = {0,139,0};
		colName2RGBdec.put("GREEN4", tmp);
		colName2RGBhex.put("GREEN4", 0xFF008B00);}
		{int[] tmp = {240,230,140};
		colName2RGBdec.put("KHAKI", tmp);
		colName2RGBhex.put("KHAKI", 0xFFF0E68C);}
		{int[] tmp = {255,246,143};
		colName2RGBdec.put("KHAKI1", tmp);
		colName2RGBhex.put("KHAKI1", 0xFFFFF68F);}
		{int[] tmp = {238,230,133};
		colName2RGBdec.put("KHAKI2", tmp);
		colName2RGBhex.put("KHAKI2", 0xFFEEE685);}
		{int[] tmp = {205,198,115};
		colName2RGBdec.put("KHAKI3", tmp);
		colName2RGBhex.put("KHAKI3", 0xFFCDC673);}
		{int[] tmp = {139,134,78};
		colName2RGBdec.put("KHAKI4", tmp);
		colName2RGBhex.put("KHAKI4", 0xFF8B864E);}
		{int[] tmp = {255,140,0};
		colName2RGBdec.put("DARKORANGE", tmp);
		colName2RGBhex.put("DARKORANGE", 0xFFFF8C00);}
		{int[] tmp = {255,127,0};
		colName2RGBdec.put("DARKORANGE1", tmp);
		colName2RGBhex.put("DARKORANGE1", 0xFFFF7F00);}
		{int[] tmp = {238,118,0};
		colName2RGBdec.put("DARKORANGE2", tmp);
		colName2RGBhex.put("DARKORANGE2", 0xFFEE7600);}
		{int[] tmp = {205,102,0};
		colName2RGBdec.put("DARKORANGE3", tmp);
		colName2RGBhex.put("DARKORANGE3", 0xFFCD6600);}
		{int[] tmp = {139,69,0};
		colName2RGBdec.put("DARKORANGE4", tmp);
		colName2RGBhex.put("DARKORANGE4", 0xFF8B4500);}
		{int[] tmp = {233,150,122};
		colName2RGBdec.put("DARKSALMON", tmp);
		colName2RGBhex.put("DARKSALMON", 0xFFE9967A);}
		{int[] tmp = {240,128,128};
		colName2RGBdec.put("LIGHTCORAL", tmp);
		colName2RGBhex.put("LIGHTCORAL", 0xFFF08080);}
		{int[] tmp = {255,160,122};
		colName2RGBdec.put("LIGHTSALMON", tmp);
		colName2RGBhex.put("LIGHTSALMON", 0xFFFFA07A);}
		{int[] tmp = {255,160,122};
		colName2RGBdec.put("LIGHTSALMON1", tmp);
		colName2RGBhex.put("LIGHTSALMON1", 0xFFFFA07A);}
		{int[] tmp = {238,149,114};
		colName2RGBdec.put("LIGHTSALMON2", tmp);
		colName2RGBhex.put("LIGHTSALMON2", 0xFFEE9572);}
		{int[] tmp = {205,129,98};
		colName2RGBdec.put("LIGHTSALMON3", tmp);
		colName2RGBhex.put("LIGHTSALMON3", 0xFFCD8162);}
		{int[] tmp = {139,87,66};
		colName2RGBdec.put("LIGHTSALMON4", tmp);
		colName2RGBhex.put("LIGHTSALMON4", 0xFF8B5742);}
		{int[] tmp = {255,218,185};
		colName2RGBdec.put("PEACHPUFF", tmp);
		colName2RGBhex.put("PEACHPUFF", 0xFFFFDAB9);}
		{int[] tmp = {255,218,185};
		colName2RGBdec.put("PEACHPUFF1", tmp);
		colName2RGBhex.put("PEACHPUFF1", 0xFFFFDAB9);}
		{int[] tmp = {238,203,173};
		colName2RGBdec.put("PEACHPUFF2", tmp);
		colName2RGBhex.put("PEACHPUFF2", 0xFFEECBAD);}
		{int[] tmp = {205,175,149};
		colName2RGBdec.put("PEACHPUFF3", tmp);
		colName2RGBhex.put("PEACHPUFF3", 0xFFCDAF95);}
		{int[] tmp = {139,119,101};
		colName2RGBdec.put("PEACHPUFF4", tmp);
		colName2RGBhex.put("PEACHPUFF4", 0xFF8B7765);}
		{int[] tmp = {255,228,196};
		colName2RGBdec.put("BISQUE", tmp);
		colName2RGBhex.put("BISQUE", 0xFFFFE4C4);}
		{int[] tmp = {255,228,196};
		colName2RGBdec.put("BISQUE1", tmp);
		colName2RGBhex.put("BISQUE1", 0xFFFFE4C4);}
		{int[] tmp = {238,213,183};
		colName2RGBdec.put("BISQUE2", tmp);
		colName2RGBhex.put("BISQUE2", 0xFFEED5B7);}
		{int[] tmp = {205,183,158};
		colName2RGBdec.put("BISQUE3", tmp);
		colName2RGBhex.put("BISQUE3", 0xFFCDB79E);}
		{int[] tmp = {139,125,107};
		colName2RGBdec.put("BISQUE4", tmp);
		colName2RGBhex.put("BISQUE4", 0xFF8B7D6B);}
		{int[] tmp = {255,127,80};
		colName2RGBdec.put("CORAL", tmp);
		colName2RGBhex.put("CORAL", 0xFFFF7F50);}
		{int[] tmp = {255,114,86};
		colName2RGBdec.put("CORAL1", tmp);
		colName2RGBhex.put("CORAL1", 0xFFFF7256);}
		{int[] tmp = {238,106,80};
		colName2RGBdec.put("CORAL2", tmp);
		colName2RGBhex.put("CORAL2", 0xFFEE6A50);}
		{int[] tmp = {205,91,69};
		colName2RGBdec.put("CORAL3", tmp);
		colName2RGBhex.put("CORAL3", 0xFFCD5B45);}
		{int[] tmp = {139,62,47};
		colName2RGBdec.put("CORAL4", tmp);
		colName2RGBhex.put("CORAL4", 0xFF8B3E2F);}
		{int[] tmp = {240,255,240};
		colName2RGBdec.put("HONEYDEW", tmp);
		colName2RGBhex.put("HONEYDEW", 0xFFF0FFF0);}
		{int[] tmp = {240,255,240};
		colName2RGBdec.put("HONEYDEW1", tmp);
		colName2RGBhex.put("HONEYDEW1", 0xFFF0FFF0);}
		{int[] tmp = {224,238,224};
		colName2RGBdec.put("HONEYDEW2", tmp);
		colName2RGBhex.put("HONEYDEW2", 0xFFE0EEE0);}
		{int[] tmp = {193,205,193};
		colName2RGBdec.put("HONEYDEW3", tmp);
		colName2RGBhex.put("HONEYDEW3", 0xFFC1CDC1);}
		{int[] tmp = {131,139,131};
		colName2RGBdec.put("HONEYDEW4", tmp);
		colName2RGBhex.put("HONEYDEW4", 0xFF838B83);}
		{int[] tmp = {255,165,0};
		colName2RGBdec.put("ORANGE", tmp);
		colName2RGBhex.put("ORANGE", 0xFFFFA500);}
		{int[] tmp = {255,165,0};
		colName2RGBdec.put("ORANGE1", tmp);
		colName2RGBhex.put("ORANGE1", 0xFFFFA500);}
		{int[] tmp = {238,154,0};
		colName2RGBdec.put("ORANGE2", tmp);
		colName2RGBhex.put("ORANGE2", 0xFFEE9A00);}
		{int[] tmp = {205,133,0};
		colName2RGBdec.put("ORANGE3", tmp);
		colName2RGBhex.put("ORANGE3", 0xFFCD8500);}
		{int[] tmp = {139,90,0};
		colName2RGBdec.put("ORANGE4", tmp);
		colName2RGBhex.put("ORANGE4", 0xFF8B5A00);}
		{int[] tmp = {250,128,114};
		colName2RGBdec.put("SALMON", tmp);
		colName2RGBhex.put("SALMON", 0xFFFA8072);}
		{int[] tmp = {255,140,105};
		colName2RGBdec.put("SALMON1", tmp);
		colName2RGBhex.put("SALMON1", 0xFFFF8C69);}
		{int[] tmp = {238,130,98};
		colName2RGBdec.put("SALMON2", tmp);
		colName2RGBhex.put("SALMON2", 0xFFEE8262);}
		{int[] tmp = {205,112,84};
		colName2RGBdec.put("SALMON3", tmp);
		colName2RGBhex.put("SALMON3", 0xFFCD7054);}
		{int[] tmp = {139,76,57};
		colName2RGBdec.put("SALMON4", tmp);
		colName2RGBhex.put("SALMON4", 0xFF8B4C39);}
		{int[] tmp = {160,82,45};
		colName2RGBdec.put("SIENNA", tmp);
		colName2RGBhex.put("SIENNA", 0xFFA0522D);}
		{int[] tmp = {255,130,71};
		colName2RGBdec.put("SIENNA1", tmp);
		colName2RGBhex.put("SIENNA1", 0xFFFF8247);}
		{int[] tmp = {238,121,66};
		colName2RGBdec.put("SIENNA2", tmp);
		colName2RGBhex.put("SIENNA2", 0xFFEE7942);}
		{int[] tmp = {205,104,57};
		colName2RGBdec.put("SIENNA3", tmp);
		colName2RGBhex.put("SIENNA3", 0xFFCD6839);}
		{int[] tmp = {139,71,38};
		colName2RGBdec.put("SIENNA4", tmp);
		colName2RGBhex.put("SIENNA4", 0xFF8B4726);}
		{int[] tmp = {255,20,147};
		colName2RGBdec.put("DEEPPINK", tmp);
		colName2RGBhex.put("DEEPPINK", 0xFFFF1493);}
		{int[] tmp = {255,20,147};
		colName2RGBdec.put("DEEPPINK1", tmp);
		colName2RGBhex.put("DEEPPINK1", 0xFFFF1493);}
		{int[] tmp = {238,18,137};
		colName2RGBdec.put("DEEPPINK2", tmp);
		colName2RGBhex.put("DEEPPINK2", 0xFFEE1289);}
		{int[] tmp = {205,16,118};
		colName2RGBdec.put("DEEPPINK3", tmp);
		colName2RGBhex.put("DEEPPINK3", 0xFFCD1076);}
		{int[] tmp = {139,10,80};
		colName2RGBdec.put("DEEPPINK4", tmp);
		colName2RGBhex.put("DEEPPINK4", 0xFF8B0A50);}
		{int[] tmp = {255,105,180};
		colName2RGBdec.put("HOTPINK", tmp);
		colName2RGBhex.put("HOTPINK", 0xFFFF69B4);}
		{int[] tmp = {255,110,180};
		colName2RGBdec.put("HOTPINK1", tmp);
		colName2RGBhex.put("HOTPINK1", 0xFFFF6EB4);}
		{int[] tmp = {238,106,167};
		colName2RGBdec.put("HOTPINK2", tmp);
		colName2RGBhex.put("HOTPINK2", 0xFFEE6AA7);}
		{int[] tmp = {205,96,144};
		colName2RGBdec.put("HOTPINK3", tmp);
		colName2RGBhex.put("HOTPINK3", 0xFFCD6090);}
		{int[] tmp = {139,58,98};
		colName2RGBdec.put("HOTPINK4", tmp);
		colName2RGBhex.put("HOTPINK4", 0xFF8B3A62);}
		{int[] tmp = {205,92,92};
		colName2RGBdec.put("INDIANRED", tmp);
		colName2RGBhex.put("INDIANRED", 0xFFCD5C5C);}
		{int[] tmp = {255,106,106};
		colName2RGBdec.put("INDIANRED1", tmp);
		colName2RGBhex.put("INDIANRED1", 0xFFFF6A6A);}
		{int[] tmp = {238,99,99};
		colName2RGBdec.put("INDIANRED2", tmp);
		colName2RGBhex.put("INDIANRED2", 0xFFEE6363);}
		{int[] tmp = {205,85,85};
		colName2RGBdec.put("INDIANRED3", tmp);
		colName2RGBhex.put("INDIANRED3", 0xFFCD5555);}
		{int[] tmp = {139,58,58};
		colName2RGBdec.put("INDIANRED4", tmp);
		colName2RGBhex.put("INDIANRED4", 0xFF8B3A3A);}
		{int[] tmp = {255,182,193};
		colName2RGBdec.put("LIGHTPINK", tmp);
		colName2RGBhex.put("LIGHTPINK", 0xFFFFB6C1);}
		{int[] tmp = {255,174,185};
		colName2RGBdec.put("LIGHTPINK1", tmp);
		colName2RGBhex.put("LIGHTPINK1", 0xFFFFAEB9);}
		{int[] tmp = {238,162,173};
		colName2RGBdec.put("LIGHTPINK2", tmp);
		colName2RGBhex.put("LIGHTPINK2", 0xFFEEA2AD);}
		{int[] tmp = {205,140,149};
		colName2RGBdec.put("LIGHTPINK3", tmp);
		colName2RGBhex.put("LIGHTPINK3", 0xFFCD8C95);}
		{int[] tmp = {139,95,101};
		colName2RGBdec.put("LIGHTPINK4", tmp);
		colName2RGBhex.put("LIGHTPINK4", 0xFF8B5F65);}
		{int[] tmp = {199,21,133};
		colName2RGBdec.put("MEDIUMVIOLETRED", tmp);
		colName2RGBhex.put("MEDIUMVIOLETRED", 0xFFC71585);}
		{int[] tmp = {255,228,225};
		colName2RGBdec.put("MISTYROSE", tmp);
		colName2RGBhex.put("MISTYROSE", 0xFFFFE4E1);}
		{int[] tmp = {255,228,225};
		colName2RGBdec.put("MISTYROSE1", tmp);
		colName2RGBhex.put("MISTYROSE1", 0xFFFFE4E1);}
		{int[] tmp = {238,213,210};
		colName2RGBdec.put("MISTYROSE2", tmp);
		colName2RGBhex.put("MISTYROSE2", 0xFFEED5D2);}
		{int[] tmp = {205,183,181};
		colName2RGBdec.put("MISTYROSE3", tmp);
		colName2RGBhex.put("MISTYROSE3", 0xFFCDB7B5);}
		{int[] tmp = {139,125,123};
		colName2RGBdec.put("MISTYROSE4", tmp);
		colName2RGBhex.put("MISTYROSE4", 0xFF8B7D7B);}
		{int[] tmp = {255,69,0};
		colName2RGBdec.put("ORANGERED", tmp);
		colName2RGBhex.put("ORANGERED", 0xFFFF4500);}
		{int[] tmp = {255,69,0};
		colName2RGBdec.put("ORANGERED1", tmp);
		colName2RGBhex.put("ORANGERED1", 0xFFFF4500);}
		{int[] tmp = {238,64,0};
		colName2RGBdec.put("ORANGERED2", tmp);
		colName2RGBhex.put("ORANGERED2", 0xFFEE4000);}
		{int[] tmp = {205,55,0};
		colName2RGBdec.put("ORANGERED3", tmp);
		colName2RGBhex.put("ORANGERED3", 0xFFCD3700);}
		{int[] tmp = {139,37,0};
		colName2RGBdec.put("ORANGERED4", tmp);
		colName2RGBhex.put("ORANGERED4", 0xFF8B2500);}
		{int[] tmp = {219,112,147};
		colName2RGBdec.put("PALEVIOLETRED", tmp);
		colName2RGBhex.put("PALEVIOLETRED", 0xFFDB7093);}
		{int[] tmp = {255,130,171};
		colName2RGBdec.put("PALEVIOLETRED1", tmp);
		colName2RGBhex.put("PALEVIOLETRED1", 0xFFFF82AB);}
		{int[] tmp = {238,121,159};
		colName2RGBdec.put("PALEVIOLETRED2", tmp);
		colName2RGBhex.put("PALEVIOLETRED2", 0xFFEE799F);}
		{int[] tmp = {205,104,137};
		colName2RGBdec.put("PALEVIOLETRED3", tmp);
		colName2RGBhex.put("PALEVIOLETRED3", 0xFFCD6889);}
		{int[] tmp = {139,71,93};
		colName2RGBdec.put("PALEVIOLETRED4", tmp);
		colName2RGBhex.put("PALEVIOLETRED4", 0xFF8B475D);}
		{int[] tmp = {208,32,144};
		colName2RGBdec.put("VIOLETRED", tmp);
		colName2RGBhex.put("VIOLETRED", 0xFFD02090);}
		{int[] tmp = {255,62,150};
		colName2RGBdec.put("VIOLETRED1", tmp);
		colName2RGBhex.put("VIOLETRED1", 0xFFFF3E96);}
		{int[] tmp = {238,58,140};
		colName2RGBdec.put("VIOLETRED2", tmp);
		colName2RGBhex.put("VIOLETRED2", 0xFFEE3A8C);}
		{int[] tmp = {205,50,120};
		colName2RGBdec.put("VIOLETRED3", tmp);
		colName2RGBhex.put("VIOLETRED3", 0xFFCD3278);}
		{int[] tmp = {139,34,82};
		colName2RGBdec.put("VIOLETRED4", tmp);
		colName2RGBhex.put("VIOLETRED4", 0xFF8B2252);}
		{int[] tmp = {178,34,34};
		colName2RGBdec.put("FIREBRICK", tmp);
		colName2RGBhex.put("FIREBRICK", 0xFFB22222);}
		{int[] tmp = {255,48,48};
		colName2RGBdec.put("FIREBRICK1", tmp);
		colName2RGBhex.put("FIREBRICK1", 0xFFFF3030);}
		{int[] tmp = {238,44,44};
		colName2RGBdec.put("FIREBRICK2", tmp);
		colName2RGBhex.put("FIREBRICK2", 0xFFEE2C2C);}
		{int[] tmp = {205,38,38};
		colName2RGBdec.put("FIREBRICK3", tmp);
		colName2RGBhex.put("FIREBRICK3", 0xFFCD2626);}
		{int[] tmp = {139,26,26};
		colName2RGBdec.put("FIREBRICK4", tmp);
		colName2RGBhex.put("FIREBRICK4", 0xFF8B1A1A);}
		{int[] tmp = {255,192,203};
		colName2RGBdec.put("PINK", tmp);
		colName2RGBhex.put("PINK", 0xFFFFC0CB);}
		{int[] tmp = {255,181,197};
		colName2RGBdec.put("PINK1", tmp);
		colName2RGBhex.put("PINK1", 0xFFFFB5C5);}
		{int[] tmp = {238,169,184};
		colName2RGBdec.put("PINK2", tmp);
		colName2RGBhex.put("PINK2", 0xFFEEA9B8);}
		{int[] tmp = {205,145,158};
		colName2RGBdec.put("PINK3", tmp);
		colName2RGBhex.put("PINK3", 0xFFCD919E);}
		{int[] tmp = {139,99,108};
		colName2RGBdec.put("PINK4", tmp);
		colName2RGBhex.put("PINK4", 0xFF8B636C);}
		{int[] tmp = {255,0,0};
		colName2RGBdec.put("RED", tmp);
		colName2RGBhex.put("RED", 0xFFFF0000);}
		{int[] tmp = {255,0,0};
		colName2RGBdec.put("RED1", tmp);
		colName2RGBhex.put("RED1", 0xFFFF0000);}
		{int[] tmp = {238,0,0};
		colName2RGBdec.put("RED2", tmp);
		colName2RGBhex.put("RED2", 0xFFEE0000);}
		{int[] tmp = {205,0,0};
		colName2RGBdec.put("RED3", tmp);
		colName2RGBhex.put("RED3", 0xFFCD0000);}
		{int[] tmp = {139,0,0};
		colName2RGBdec.put("RED4", tmp);
		colName2RGBhex.put("RED4", 0xFF8B0000);}
		{int[] tmp = {255,99,71};
		colName2RGBdec.put("TOMATO", tmp);
		colName2RGBhex.put("TOMATO", 0xFFFF6347);}
		{int[] tmp = {255,99,71};
		colName2RGBdec.put("TOMATO1", tmp);
		colName2RGBhex.put("TOMATO1", 0xFFFF6347);}
		{int[] tmp = {238,92,66};
		colName2RGBdec.put("TOMATO2", tmp);
		colName2RGBhex.put("TOMATO2", 0xFFEE5C42);}
		{int[] tmp = {205,79,57};
		colName2RGBdec.put("TOMATO3", tmp);
		colName2RGBhex.put("TOMATO3", 0xFFCD4F39);}
		{int[] tmp = {139,54,38};
		colName2RGBdec.put("TOMATO4", tmp);
		colName2RGBhex.put("TOMATO4", 0xFF8B3626);}
		{int[] tmp = {153,50,204};
		colName2RGBdec.put("DARKORCHID", tmp);
		colName2RGBhex.put("DARKORCHID", 0xFF9932CC);}
		{int[] tmp = {191,62,255};
		colName2RGBdec.put("DARKORCHID1", tmp);
		colName2RGBhex.put("DARKORCHID1", 0xFFBF3EFF);}
		{int[] tmp = {178,58,238};
		colName2RGBdec.put("DARKORCHID2", tmp);
		colName2RGBhex.put("DARKORCHID2", 0xFFB23AEE);}
		{int[] tmp = {154,50,205};
		colName2RGBdec.put("DARKORCHID3", tmp);
		colName2RGBhex.put("DARKORCHID3", 0xFF9A32CD);}
		{int[] tmp = {104,34,139};
		colName2RGBdec.put("DARKORCHID4", tmp);
		colName2RGBhex.put("DARKORCHID4", 0xFF68228B);}
		{int[] tmp = {148,0,211};
		colName2RGBdec.put("DARKVIOLET", tmp);
		colName2RGBhex.put("DARKVIOLET", 0xFF9400D3);}
		{int[] tmp = {255,240,245};
		colName2RGBdec.put("LAVENDERBLUSH", tmp);
		colName2RGBhex.put("LAVENDERBLUSH", 0xFFFFF0F5);}
		{int[] tmp = {255,240,245};
		colName2RGBdec.put("LAVENDERBLUSH1", tmp);
		colName2RGBhex.put("LAVENDERBLUSH1", 0xFFFFF0F5);}
		{int[] tmp = {238,224,229};
		colName2RGBdec.put("LAVENDERBLUSH2", tmp);
		colName2RGBhex.put("LAVENDERBLUSH2", 0xFFEEE0E5);}
		{int[] tmp = {205,193,197};
		colName2RGBdec.put("LAVENDERBLUSH3", tmp);
		colName2RGBhex.put("LAVENDERBLUSH3", 0xFFCDC1C5);}
		{int[] tmp = {139,131,134};
		colName2RGBdec.put("LAVENDERBLUSH4", tmp);
		colName2RGBhex.put("LAVENDERBLUSH4", 0xFF8B8386);}
		{int[] tmp = {186,85,211};
		colName2RGBdec.put("MEDIUMORCHID", tmp);
		colName2RGBhex.put("MEDIUMORCHID", 0xFFBA55D3);}
		{int[] tmp = {224,102,255};
		colName2RGBdec.put("MEDIUMORCHID1", tmp);
		colName2RGBhex.put("MEDIUMORCHID1", 0xFFE066FF);}
		{int[] tmp = {209,95,238};
		colName2RGBdec.put("MEDIUMORCHID2", tmp);
		colName2RGBhex.put("MEDIUMORCHID2", 0xFFD15FEE);}
		{int[] tmp = {180,82,205};
		colName2RGBdec.put("MEDIUMORCHID3", tmp);
		colName2RGBhex.put("MEDIUMORCHID3", 0xFFB452CD);}
		{int[] tmp = {122,55,139};
		colName2RGBdec.put("MEDIUMORCHID4", tmp);
		colName2RGBhex.put("MEDIUMORCHID4", 0xFF7A378B);}
		{int[] tmp = {147,112,219};
		colName2RGBdec.put("MEDIUMPURPLE", tmp);
		colName2RGBhex.put("MEDIUMPURPLE", 0xFF9370DB);}
		{int[] tmp = {171,130,255};
		colName2RGBdec.put("MEDIUMPURPLE1", tmp);
		colName2RGBhex.put("MEDIUMPURPLE1", 0xFFAB82FF);}
		{int[] tmp = {159,121,238};
		colName2RGBdec.put("MEDIUMPURPLE2", tmp);
		colName2RGBhex.put("MEDIUMPURPLE2", 0xFF9F79EE);}
		{int[] tmp = {137,104,205};
		colName2RGBdec.put("MEDIUMPURPLE3", tmp);
		colName2RGBhex.put("MEDIUMPURPLE3", 0xFF8968CD);}
		{int[] tmp = {93,71,139};
		colName2RGBdec.put("MEDIUMPURPLE4", tmp);
		colName2RGBhex.put("MEDIUMPURPLE4", 0xFF5D478B);}
		{int[] tmp = {230,230,250};
		colName2RGBdec.put("LAVENDER", tmp);
		colName2RGBhex.put("LAVENDER", 0xFFE6E6FA);}
		{int[] tmp = {255,0,255};
		colName2RGBdec.put("MAGENTA", tmp);
		colName2RGBhex.put("MAGENTA", 0xFFFF00FF);}
		{int[] tmp = {255,0,255};
		colName2RGBdec.put("MAGENTA1", tmp);
		colName2RGBhex.put("MAGENTA1", 0xFFFF00FF);}
		{int[] tmp = {238,0,238};
		colName2RGBdec.put("MAGENTA2", tmp);
		colName2RGBhex.put("MAGENTA2", 0xFFEE00EE);}
		{int[] tmp = {205,0,205};
		colName2RGBdec.put("MAGENTA3", tmp);
		colName2RGBhex.put("MAGENTA3", 0xFFCD00CD);}
		{int[] tmp = {139,0,139};
		colName2RGBdec.put("MAGENTA4", tmp);
		colName2RGBhex.put("MAGENTA4", 0xFF8B008B);}
		{int[] tmp = {176,48,96};
		colName2RGBdec.put("MAROON", tmp);
		colName2RGBhex.put("MAROON", 0xFFB03060);}
		{int[] tmp = {255,52,179};
		colName2RGBdec.put("MAROON1", tmp);
		colName2RGBhex.put("MAROON1", 0xFFFF34B3);}
		{int[] tmp = {238,48,167};
		colName2RGBdec.put("MAROON2", tmp);
		colName2RGBhex.put("MAROON2", 0xFFEE30A7);}
		{int[] tmp = {205,41,144};
		colName2RGBdec.put("MAROON3", tmp);
		colName2RGBhex.put("MAROON3", 0xFFCD2990);}
		{int[] tmp = {139,28,98};
		colName2RGBdec.put("MAROON4", tmp);
		colName2RGBhex.put("MAROON4", 0xFF8B1C62);}
		{int[] tmp = {218,112,214};
		colName2RGBdec.put("ORCHID", tmp);
		colName2RGBhex.put("ORCHID", 0xFFDA70D6);}
		{int[] tmp = {255,131,250};
		colName2RGBdec.put("ORCHID1", tmp);
		colName2RGBhex.put("ORCHID1", 0xFFFF83FA);}
		{int[] tmp = {238,122,233};
		colName2RGBdec.put("ORCHID2", tmp);
		colName2RGBhex.put("ORCHID2", 0xFFEE7AE9);}
		{int[] tmp = {205,105,201};
		colName2RGBdec.put("ORCHID3", tmp);
		colName2RGBhex.put("ORCHID3", 0xFFCD69C9);}
		{int[] tmp = {139,71,137};
		colName2RGBdec.put("ORCHID4", tmp);
		colName2RGBhex.put("ORCHID4", 0xFF8B4789);}
		{int[] tmp = {221,160,221};
		colName2RGBdec.put("PLUM", tmp);
		colName2RGBhex.put("PLUM", 0xFFDDA0DD);}
		{int[] tmp = {255,187,255};
		colName2RGBdec.put("PLUM1", tmp);
		colName2RGBhex.put("PLUM1", 0xFFFFBBFF);}
		{int[] tmp = {238,174,238};
		colName2RGBdec.put("PLUM2", tmp);
		colName2RGBhex.put("PLUM2", 0xFFEEAEEE);}
		{int[] tmp = {205,150,205};
		colName2RGBdec.put("PLUM3", tmp);
		colName2RGBhex.put("PLUM3", 0xFFCD96CD);}
		{int[] tmp = {139,102,139};
		colName2RGBdec.put("PLUM4", tmp);
		colName2RGBhex.put("PLUM4", 0xFF8B668B);}
		{int[] tmp = {160,32,240};
		colName2RGBdec.put("PURPLE", tmp);
		colName2RGBhex.put("PURPLE", 0xFFA020F0);}
		{int[] tmp = {155,48,255};
		colName2RGBdec.put("PURPLE1", tmp);
		colName2RGBhex.put("PURPLE1", 0xFF9B30FF);}
		{int[] tmp = {145,44,238};
		colName2RGBdec.put("PURPLE2", tmp);
		colName2RGBhex.put("PURPLE2", 0xFF912CEE);}
		{int[] tmp = {125,38,205};
		colName2RGBdec.put("PURPLE3", tmp);
		colName2RGBhex.put("PURPLE3", 0xFF7D26CD);}
		{int[] tmp = {85,26,139};
		colName2RGBdec.put("PURPLE4", tmp);
		colName2RGBhex.put("PURPLE4", 0xFF551A8B);}
		{int[] tmp = {216,191,216};
		colName2RGBdec.put("THISTLE", tmp);
		colName2RGBhex.put("THISTLE", 0xFFD8BFD8);}
		{int[] tmp = {255,225,255};
		colName2RGBdec.put("THISTLE1", tmp);
		colName2RGBhex.put("THISTLE1", 0xFFFFE1FF);}
		{int[] tmp = {238,210,238};
		colName2RGBdec.put("THISTLE2", tmp);
		colName2RGBhex.put("THISTLE2", 0xFFEED2EE);}
		{int[] tmp = {205,181,205};
		colName2RGBdec.put("THISTLE3", tmp);
		colName2RGBhex.put("THISTLE3", 0xFFCDB5CD);}
		{int[] tmp = {139,123,139};
		colName2RGBdec.put("THISTLE4", tmp);
		colName2RGBhex.put("THISTLE4", 0xFF8B7B8B);}
		{int[] tmp = {238,130,238};
		colName2RGBdec.put("VIOLET", tmp);
		colName2RGBhex.put("VIOLET", 0xFFEE82EE);}
		{int[] tmp = {250,235,215};
		colName2RGBdec.put("ANTIQUEWHITE", tmp);
		colName2RGBhex.put("ANTIQUEWHITE", 0xFFFAEBD7);}
		{int[] tmp = {255,239,219};
		colName2RGBdec.put("ANTIQUEWHITE1", tmp);
		colName2RGBhex.put("ANTIQUEWHITE1", 0xFFFFEFDB);}
		{int[] tmp = {238,223,204};
		colName2RGBdec.put("ANTIQUEWHITE2", tmp);
		colName2RGBhex.put("ANTIQUEWHITE2", 0xFFEEDFCC);}
		{int[] tmp = {205,192,176};
		colName2RGBdec.put("ANTIQUEWHITE3", tmp);
		colName2RGBhex.put("ANTIQUEWHITE3", 0xFFCDC0B0);}
		{int[] tmp = {139,131,120};
		colName2RGBdec.put("ANTIQUEWHITE4", tmp);
		colName2RGBhex.put("ANTIQUEWHITE4", 0xFF8B8378);}
		{int[] tmp = {255,250,240};
		colName2RGBdec.put("FLORALWHITE", tmp);
		colName2RGBhex.put("FLORALWHITE", 0xFFFFFAF0);}
		{int[] tmp = {248,248,255};
		colName2RGBdec.put("GHOSTWHITE", tmp);
		colName2RGBhex.put("GHOSTWHITE", 0xFFF8F8FF);}
		{int[] tmp = {255,222,173};
		colName2RGBdec.put("NAVAJOWHITE", tmp);
		colName2RGBhex.put("NAVAJOWHITE", 0xFFFFDEAD);}
		{int[] tmp = {255,222,173};
		colName2RGBdec.put("NAVAJOWHITE1", tmp);
		colName2RGBhex.put("NAVAJOWHITE1", 0xFFFFDEAD);}
		{int[] tmp = {238,207,161};
		colName2RGBdec.put("NAVAJOWHITE2", tmp);
		colName2RGBhex.put("NAVAJOWHITE2", 0xFFEECFA1);}
		{int[] tmp = {205,179,139};
		colName2RGBdec.put("NAVAJOWHITE3", tmp);
		colName2RGBhex.put("NAVAJOWHITE3", 0xFFCDB38B);}
		{int[] tmp = {139,121,94};
		colName2RGBdec.put("NAVAJOWHITE4", tmp);
		colName2RGBhex.put("NAVAJOWHITE4", 0xFF8B795E);}
		{int[] tmp = {253,245,230};
		colName2RGBdec.put("OLDLACE", tmp);
		colName2RGBhex.put("OLDLACE", 0xFFFDF5E6);}
		{int[] tmp = {245,245,245};
		colName2RGBdec.put("WHITESMOKE", tmp);
		colName2RGBhex.put("WHITESMOKE", 0xFFF5F5F5);}
		{int[] tmp = {220,220,220};
		colName2RGBdec.put("GAINSBORO", tmp);
		colName2RGBhex.put("GAINSBORO", 0xFFDCDCDC);}
		{int[] tmp = {255,255,240};
		colName2RGBdec.put("IVORY", tmp);
		colName2RGBhex.put("IVORY", 0xFFFFFFF0);}
		{int[] tmp = {255,255,240};
		colName2RGBdec.put("IVORY1", tmp);
		colName2RGBhex.put("IVORY1", 0xFFFFFFF0);}
		{int[] tmp = {238,238,224};
		colName2RGBdec.put("IVORY2", tmp);
		colName2RGBhex.put("IVORY2", 0xFFEEEEE0);}
		{int[] tmp = {205,205,193};
		colName2RGBdec.put("IVORY3", tmp);
		colName2RGBhex.put("IVORY3", 0xFFCDCDC1);}
		{int[] tmp = {139,139,131};
		colName2RGBdec.put("IVORY4", tmp);
		colName2RGBhex.put("IVORY4", 0xFF8B8B83);}
		{int[] tmp = {250,240,230};
		colName2RGBdec.put("LINEN", tmp);
		colName2RGBhex.put("LINEN", 0xFFFAF0E6);}
		{int[] tmp = {255,245,238};
		colName2RGBdec.put("SEASHELL", tmp);
		colName2RGBhex.put("SEASHELL", 0xFFFFF5EE);}
		{int[] tmp = {255,245,238};
		colName2RGBdec.put("SEASHELL1", tmp);
		colName2RGBhex.put("SEASHELL1", 0xFFFFF5EE);}
		{int[] tmp = {238,229,222};
		colName2RGBdec.put("SEASHELL2", tmp);
		colName2RGBhex.put("SEASHELL2", 0xFFEEE5DE);}
		{int[] tmp = {205,197,191};
		colName2RGBdec.put("SEASHELL3", tmp);
		colName2RGBhex.put("SEASHELL3", 0xFFCDC5BF);}
		{int[] tmp = {139,134,130};
		colName2RGBdec.put("SEASHELL4", tmp);
		colName2RGBhex.put("SEASHELL4", 0xFF8B8682);}
		{int[] tmp = {255,250,250};
		colName2RGBdec.put("SNOW", tmp);
		colName2RGBhex.put("SNOW", 0xFFFFFAFA);}
		{int[] tmp = {255,250,250};
		colName2RGBdec.put("SNOW1", tmp);
		colName2RGBhex.put("SNOW1", 0xFFFFFAFA);}
		{int[] tmp = {238,233,233};
		colName2RGBdec.put("SNOW2", tmp);
		colName2RGBhex.put("SNOW2", 0xFFEEE9E9);}
		{int[] tmp = {205,201,201};
		colName2RGBdec.put("SNOW3", tmp);
		colName2RGBhex.put("SNOW3", 0xFFCDC9C9);}
		{int[] tmp = {139,137,137};
		colName2RGBdec.put("SNOW4", tmp);
		colName2RGBhex.put("SNOW4", 0xFF8B8989);}
		{int[] tmp = {245,222,179};
		colName2RGBdec.put("WHEAT", tmp);
		colName2RGBhex.put("WHEAT", 0xFFF5DEB3);}
		{int[] tmp = {255,231,186};
		colName2RGBdec.put("WHEAT1", tmp);
		colName2RGBhex.put("WHEAT1", 0xFFFFE7BA);}
		{int[] tmp = {238,216,174};
		colName2RGBdec.put("WHEAT2", tmp);
		colName2RGBhex.put("WHEAT2", 0xFFEED8AE);}
		{int[] tmp = {205,186,150};
		colName2RGBdec.put("WHEAT3", tmp);
		colName2RGBhex.put("WHEAT3", 0xFFCDBA96);}
		{int[] tmp = {139,126,102};
		colName2RGBdec.put("WHEAT4", tmp);
		colName2RGBhex.put("WHEAT4", 0xFF8B7E66);}
		{int[] tmp = {255,235,205};
		colName2RGBdec.put("BLANCHEDALMOND", tmp);
		colName2RGBhex.put("BLANCHEDALMOND", 0xFFFFEBCD);}
		{int[] tmp = {184,134,11};
		colName2RGBdec.put("DARKGOLDENROD", tmp);
		colName2RGBhex.put("DARKGOLDENROD", 0xFFB8860B);}
		{int[] tmp = {255,185,15};
		colName2RGBdec.put("DARKGOLDENROD1", tmp);
		colName2RGBhex.put("DARKGOLDENROD1", 0xFFFFB90F);}
		{int[] tmp = {238,173,14};
		colName2RGBdec.put("DARKGOLDENROD2", tmp);
		colName2RGBhex.put("DARKGOLDENROD2", 0xFFEEAD0E);}
		{int[] tmp = {205,149,12};
		colName2RGBdec.put("DARKGOLDENROD3", tmp);
		colName2RGBhex.put("DARKGOLDENROD3", 0xFFCD950C);}
		{int[] tmp = {139,101,8};
		colName2RGBdec.put("DARKGOLDENROD4", tmp);
		colName2RGBhex.put("DARKGOLDENROD4", 0xFF8B6508);}
		{int[] tmp = {255,250,205};
		colName2RGBdec.put("LEMONCHIFFON", tmp);
		colName2RGBhex.put("LEMONCHIFFON", 0xFFFFFACD);}
		{int[] tmp = {255,250,205};
		colName2RGBdec.put("LEMONCHIFFON1", tmp);
		colName2RGBhex.put("LEMONCHIFFON1", 0xFFFFFACD);}
		{int[] tmp = {238,233,191};
		colName2RGBdec.put("LEMONCHIFFON2", tmp);
		colName2RGBhex.put("LEMONCHIFFON2", 0xFFEEE9BF);}
		{int[] tmp = {205,201,165};
		colName2RGBdec.put("LEMONCHIFFON3", tmp);
		colName2RGBhex.put("LEMONCHIFFON3", 0xFFCDC9A5);}
		{int[] tmp = {139,137,112};
		colName2RGBdec.put("LEMONCHIFFON4", tmp);
		colName2RGBhex.put("LEMONCHIFFON4", 0xFF8B8970);}
		{int[] tmp = {238,221,130};
		colName2RGBdec.put("LIGHTGOLDENROD", tmp);
		colName2RGBhex.put("LIGHTGOLDENROD", 0xFFEEDD82);}
		{int[] tmp = {255,236,139};
		colName2RGBdec.put("LIGHTGOLDENROD1", tmp);
		colName2RGBhex.put("LIGHTGOLDENROD1", 0xFFFFEC8B);}
		{int[] tmp = {238,220,130};
		colName2RGBdec.put("LIGHTGOLDENROD2", tmp);
		colName2RGBhex.put("LIGHTGOLDENROD2", 0xFFEEDC82);}
		{int[] tmp = {205,190,112};
		colName2RGBdec.put("LIGHTGOLDENROD3", tmp);
		colName2RGBhex.put("LIGHTGOLDENROD3", 0xFFCDBE70);}
		{int[] tmp = {139,129,76};
		colName2RGBdec.put("LIGHTGOLDENROD4", tmp);
		colName2RGBhex.put("LIGHTGOLDENROD4", 0xFF8B814C);}
		{int[] tmp = {250,250,210};
		colName2RGBdec.put("LIGHTGOLDENRODYELLOW", tmp);
		colName2RGBhex.put("LIGHTGOLDENRODYELLOW", 0xFFFAFAD2);}
		{int[] tmp = {255,255,224};
		colName2RGBdec.put("LIGHTYELLOW", tmp);
		colName2RGBhex.put("LIGHTYELLOW", 0xFFFFFFE0);}
		{int[] tmp = {255,255,224};
		colName2RGBdec.put("LIGHTYELLOW1", tmp);
		colName2RGBhex.put("LIGHTYELLOW1", 0xFFFFFFE0);}
		{int[] tmp = {238,238,209};
		colName2RGBdec.put("LIGHTYELLOW2", tmp);
		colName2RGBhex.put("LIGHTYELLOW2", 0xFFEEEED1);}
		{int[] tmp = {205,205,180};
		colName2RGBdec.put("LIGHTYELLOW3", tmp);
		colName2RGBhex.put("LIGHTYELLOW3", 0xFFCDCDB4);}
		{int[] tmp = {139,139,122};
		colName2RGBdec.put("LIGHTYELLOW4", tmp);
		colName2RGBhex.put("LIGHTYELLOW4", 0xFF8B8B7A);}
		{int[] tmp = {238,232,170};
		colName2RGBdec.put("PALEGOLDENROD", tmp);
		colName2RGBhex.put("PALEGOLDENROD", 0xFFEEE8AA);}
		{int[] tmp = {255,239,213};
		colName2RGBdec.put("PAPAYAWHIP", tmp);
		colName2RGBhex.put("PAPAYAWHIP", 0xFFFFEFD5);}
		{int[] tmp = {255,248,220};
		colName2RGBdec.put("CORNSILK", tmp);
		colName2RGBhex.put("CORNSILK", 0xFFFFF8DC);}
		{int[] tmp = {255,248,220};
		colName2RGBdec.put("CORNSILK1", tmp);
		colName2RGBhex.put("CORNSILK1", 0xFFFFF8DC);}
		{int[] tmp = {238,232,205};
		colName2RGBdec.put("CORNSILK2", tmp);
		colName2RGBhex.put("CORNSILK2", 0xFFEEE8CD);}
		{int[] tmp = {205,200,177};
		colName2RGBdec.put("CORNSILK3", tmp);
		colName2RGBhex.put("CORNSILK3", 0xFFCDC8B1);}
		{int[] tmp = {139,136,120};
		colName2RGBdec.put("CORNSILK4", tmp);
		colName2RGBhex.put("CORNSILK4", 0xFF8B8878);}
		{int[] tmp = {218,165,32};
		colName2RGBdec.put("GOLDENROD", tmp);
		colName2RGBhex.put("GOLDENROD", 0xFFDAA520);}
		{int[] tmp = {255,193,37};
		colName2RGBdec.put("GOLDENROD1", tmp);
		colName2RGBhex.put("GOLDENROD1", 0xFFFFC125);}
		{int[] tmp = {238,180,34};
		colName2RGBdec.put("GOLDENROD2", tmp);
		colName2RGBhex.put("GOLDENROD2", 0xFFEEB422);}
		{int[] tmp = {205,155,29};
		colName2RGBdec.put("GOLDENROD3", tmp);
		colName2RGBhex.put("GOLDENROD3", 0xFFCD9B1D);}
		{int[] tmp = {139,105,20};
		colName2RGBdec.put("GOLDENROD4", tmp);
		colName2RGBhex.put("GOLDENROD4", 0xFF8B6914);}
		{int[] tmp = {255,228,181};
		colName2RGBdec.put("MOCCASIN", tmp);
		colName2RGBhex.put("MOCCASIN", 0xFFFFE4B5);}
		{int[] tmp = {255,255,0};
		colName2RGBdec.put("YELLOW", tmp);
		colName2RGBhex.put("YELLOW", 0xFFFFFF00);}
		{int[] tmp = {255,255,0};
		colName2RGBdec.put("YELLOW1", tmp);
		colName2RGBhex.put("YELLOW1", 0xFFFFFF00);}
		{int[] tmp = {238,238,0};
		colName2RGBdec.put("YELLOW2", tmp);
		colName2RGBhex.put("YELLOW2", 0xFFEEEE00);}
		{int[] tmp = {205,205,0};
		colName2RGBdec.put("YELLOW3", tmp);
		colName2RGBhex.put("YELLOW3", 0xFFCDCD00);}
		{int[] tmp = {139,139,0};
		colName2RGBdec.put("YELLOW4", tmp);
		colName2RGBhex.put("YELLOW4", 0xFF8B8B00);}
		{int[] tmp = {255,215,0};
		colName2RGBdec.put("GOLD", tmp);
		colName2RGBhex.put("GOLD", 0xFFFFD700);}
		{int[] tmp = {255,215,0};
		colName2RGBdec.put("GOLD1", tmp);
		colName2RGBhex.put("GOLD1", 0xFFFFD700);}
		{int[] tmp = {238,201,0};
		colName2RGBdec.put("GOLD2", tmp);
		colName2RGBhex.put("GOLD2", 0xFFEEC900);}
		{int[] tmp = {205,173,0};
		colName2RGBdec.put("GOLD3", tmp);
		colName2RGBhex.put("GOLD3", 0xFFCDAD00);}
		{int[] tmp = {139,117,0};
		colName2RGBdec.put("GOLD4", tmp);
		colName2RGBhex.put("GOLD4", 0xFF8B7500);}


	}
	
	private static String defaultColor  = "WHITE";
	
	public static boolean contains(String colName){
		return colName2RGBdec.containsKey(colName.toUpperCase());
	}
	
	public static int[] getRGBdec(String colName){
		String tmp  = colName.toUpperCase();
		if(colName2RGBdec.containsKey(tmp)){
			return colName2RGBdec.get(tmp);
		}else{
			return colName2RGBdec.get(defaultColor);
		}
	}
	public static int getRGBhex(String colName){
		String tmp  = colName.toUpperCase();
		if(colName2RGBhex.containsKey(tmp)){
			return colName2RGBhex.get(tmp);	
		}else{
			return colName2RGBhex.get(defaultColor);
		}
	}
	public static int lerpColor(int c1, int c2, int c3, float amt, int mode){
		if(amt < 0.5){
			return PApplet.lerpColor(c1, c2, amt*2, mode);
		}else{
			return PApplet.lerpColor(c2, c3, (float)(amt-0.5)*2, mode);
		}  	
	}
}
