/*
 * MIT License
 * 
 * Copyright (c) 2021 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
 *                    Nils Hoffmann  -  nils.hoffmann {at} isas.de
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the 'Software'), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:;
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHether IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

grammar FattyAcids;


/* first rule is always start rule, EOF = end of file */
lipid : lipid_eof EOF;


lipid_eof : fatty_acid;
fatty_acid: regular_fatty_acid_types | wax | CAR | ethanolamine;
wax : wax_ester SPACE fatty_acid_type | wax_ester SPACE additional_descriptions fatty_acid_type;
regular_fatty_acid_types : regular_fatty_acid | regular_fatty_acid cyclo | regular_fatty_acid CoA;

regular_fatty_acid : fatty_acid_type | additional_descriptions fatty_acid_type | additional_descriptions DASH ROB double_bond_positions RCB DASH fatty_length acid_type_double | ROB double_bond_positions RCB DASH fatty_length acid_type_double;



fatty_acid_type : acid_types | cycle acid_types | ate_type;
acid_types : fatty_length acid_type_regular | double_bond_positions DASH fatty_length acid_type_double | fatty_length acid_type_double | fatty_length acid_type_double;

fatty_length : notation_specials | notation_regular;
notation_regular : notation_last_digit | notation_last_digit notation_second_digit | notation_second_digit;
/* 1, 2, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9 */
notation_last_digit : 'un' | 'hen' | 'do' | 'di' | 'tri' | 'buta' | 'but' | 'tetra' | 'penta' | 'pent' | 'hexa' | 'hex' | 'hepta' | 'hept' | 'octa' | 'oct' | 'nona' | 'non';
/* 0, 10, 10, 20, 20, 30 */
notation_second_digit: 'deca' | 'dec' | 'cosa' | 'cos' | 'eicosa' | 'eicos' | 'triaconta' | 'triacont' | 'tetraconta'  | 'tetracont' | 'pentaconta' | 'pantacont';
/* 4, 10, 20, 21, 21, 30, 30 */
notation_specials: 'etha' | 'eth' | 'buta' | 'but' | 'butr' | 'valer' | 'propa' | 'propi' | 'propio' | 'prop' | 'eicosa' | 'eicos' | 'icosa' | 'icos' | 'heneicosa' | 'heneicos' | prosta | isoprop;
isoprop: 'isoprop';
prosta : 'prosta' | 'prost' | 'prostan';
CoA : DASH 'coa';

acid_type_regular: acid_single_type;
acid_type_double: db_num acid_single_type;
acid_single_type: 'noic acid' | 'nal' | dioic | 'noyloxy' | '-1-yl' | 'noyl' | 'nyl' | 'yl' | ol | dial | 'noate' | 'nate'; /* | 'ne' */

db_num: DASH double_bond_positions DASH db_length db_suffix | DASH double_bond_positions DASH db_suffix | db_length db_suffix | db_suffix;
db_suffix : 'e' | 'ne' | 'ene' | 'en' | 'n';
dial : 'dial';
db_length: notation_regular;
dioic : 'n' DASH functional_pos pos_separator functional_pos DASH dioic_acid | 'ne' dioic_acid;
dioic_acid : 'dioic acid';
ol : 'nol' | db_suffix DASH hydroxyl_positions DASH notation_regular 'ol' | db_suffix DASH hydroxyl_position DASH 'ol' | DASH hydroxyl_positions DASH notation_regular 'ol' | DASH hydroxyl_position DASH 'ol';
ate_type : ate | additional_descriptions ate;
ate : 'formate' | 'acetate' | 'butyrate' | 'propionate' | 'valerate' | isobut;
isobut : 'isobutyrate';


hydroxyl_positions : hydroxyl_positions pos_separator hydroxyl_positions | hydroxyl_position;
hydroxyl_position : hydroxyl_number | hydroxyl_number stereo | hydroxyl_number '\'' stereo;
hydroxyl_number : number;

additional_descriptions : additional_descriptions_m | fg_pos_summary DASH additional_descriptions_m | fg_pos_summary DASH;
additional_descriptions_m : additional_descriptions_m additional_descriptions_m | additional_description;
additional_description : functional_group | functional_group DASH | pos_neg | reduction | reduction DASH;
functional_group : multi_functional_group | single_functional_group | epoxy;
pos_neg : '(+/-)-' | '(+)-' | '(-)-';
double_bond_positions : double_bond_positions pos_separator double_bond_positions | double_bond_position;
double_bond_position : db_number | db_number cistrans | db_number '\'' | db_number '\'' cistrans | cistrans;
cistrans : 'e' | 'z' | 'r' | 's' | 'a' | 'b' | 'c';
db_number : number;
fg_pos_summary : ROB functional_positions RCB;

multi_functional_group : functional_positions DASH functional_length functional_group_type | functional_positions DASH functional_group_type;
functional_length : notation_last_digit | notation_second_digit | notation_last_digit notation_second_digit;
functional_positions : functional_positions_pure | ROB functional_positions_pure RCB;
functional_positions_pure : functional_positions pos_separator functional_positions | functional_position;
single_functional_group : functional_position DASH functional_group_type_name | functional_position functional_group_type_name | recursion_description | recursion_description DASH;
functional_group_type_name : functional_group_type | ROB functional_group_type RCB;
functional_group_type : 'hydroxy' | 'oxo' | 'bromo' | 'thio' | 'keto' | 'methyl' | 'hydroperoxy' | homo | 'fluoro' | 'chloro' | methylene | 'sulfooxy' | 'amino' | 'sulfanyl' | 'methoxy' | 'iodo' | 'cyano' | 'nitro' | 'OH' | 'thio' | 'mercapto' | 'carboxy';
epoxy : functional_position pos_separator functional_position DASH 'epoxy' | functional_position ROB functional_position RCB DASH 'epoxy';
methylene : 'methylene';

functional_position : functional_position_pure | ROB functional_position_pure RCB;
functional_position_pure : functional_pos_pr | functional_pos func_stereo | func_stereo;
functional_pos_pr : functional_pos | functional_pos '\'';
functional_pos : number;
func_stereo : stereo;
stereo : cistrans;
reduction : functional_position DASH 'nor' | functional_positions DASH functional_length 'nor';
homo : 'homo';

cycle : 'cyclo';
cyclo : '-cyclo' SOB functional_position pos_separator functional_position SCB;

recursion_description : recursion_position DASH recursion | recursion_position DASH ROB recursion RCB;
recursion : fatty_acid;
recursion_position : ROB functional_positions RCB | recursion_pos | recursion_pos stereo;
recursion_pos : number;

wax_ester : fatty_acid | ROB fatty_acid RCB | methyl;
methyl : 'methyl';
CAR : car_positions DASH SOB fatty_acid SCB '-4-(trimethylazaniumyl)butanoate';
car_positions : functional_position | ROB car_position RCB DASH functional_position;
ethanolamine : 'n-' ROB fatty_acid RCB DASH 'ethanolamine';



/* separators */
SPACE : ' ';
COLON : ':';
SEMICOLON : ';';
DASH : '-' | '‐';
UNDERSCORE : '_';
SLASH : '/';
BACKSLASH : '\\';
COMMA: ',';
ROB: '(';
RCB: ')';
SOB: '[';
SCB: ']';
pos_separator : COMMA;



number :  digit;
digit : '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | digit digit;


