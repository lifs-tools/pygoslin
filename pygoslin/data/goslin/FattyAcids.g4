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
fatty_acid: fatty_length acid_description | additional_descriptions fatty_length acid_description;
acid_description : acid_type | acid_type cyclo;

fatty_length : notation_specials | notation_last_digit | notation_last_digit notation_second_digit;
/* 1, 2, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9 */
notation_last_digit : 'un' | 'do' | 'di' | 'tri' | 'buta' | 'but' | 'tetra' | 'penta' | 'pent' | 'hexa' | 'hex' | 'hepta' | 'hept' | 'octa' | 'oct' | 'nona' | 'non';
/* 0, 10, 10, 20, 20, 30 */
notation_second_digit: 'deca' | 'dec' | 'cosa' | 'cos' | 'triaconta' | 'triacont' | 'tetraconta'  | 'tetracont';
/* 4, 10, 20, 21, 21, 30, 30 */
notation_specials: 'buta' | 'but' | 'deca' | 'dec' | 'eicosa' | 'eicos' | 'heneicosa' | 'heneicos' | 'triaconta' | 'triacont' | 'tetraconta' | 'tetracont' | prosta;
prosta : 'prosta' | 'prost' | 'prostan';

acid_type: db_num acid_single_type | acid_single_type;
acid_single_type: 'nol' | 'noic acid' | 'nal' | dioic | 'nyl acetate' | 'noyloxy';
db_num: DASH double_bond_positions DASH db_length db_suffix | db_length db_suffix | db_suffix;
db_suffix : 'e' | 'ne' | 'ene';
db_length: notation_last_digit;
dioic : DASH functional_pos pos_separator functional_pos DASH dioic_acid | dioic_acid;
dioic_acid : 'dioic acid';

additional_descriptions : additional_descriptions additional_descriptions | additional_description;
additional_description : functional_group | functional_group DASH | double_bond_positions DASH | ROB double_bond_positions RCB DASH | '(+/-)-' | reduction | reduction DASH | fahfa_description;
functional_group : multi_functional_group | single_functional_group | epoxy;
double_bond_positions : double_bond_positions pos_separator double_bond_positions | double_bond_position;
double_bond_position : db_number | db_number cistrans;
cistrans : 'e' | 'z';
db_number : number;

multi_functional_group : functional_positions DASH functional_length functional_group_type;
functional_length : notation_last_digit | notation_second_digit | notation_last_digit notation_second_digit;
functional_positions : functional_positions pos_separator functional_positions | functional_position;
single_functional_group : functional_position DASH functional_group_type | functional_position functional_group_type;
functional_group_type : 'ethyl' | 'propyl' | 'hydroxy' | 'oxo' | 'bromo' | 'thio' | 'keto' | 'methyl' | 'hydroperoxy' | 'homo' | 'fluoro';
epoxy : functional_position pos_separator functional_position DASH 'epoxy' | functional_position ROB functional_position RCB DASH 'epoxy';

functional_position : functional_pos | functional_pos stereo;
functional_pos : number;
stereo : 'r' | 's' | 'a' | 'b';
reduction : functional_position DASH 'nor' | functional_positions DASH functional_length 'nor';

cyclo : '-cyclo' SOB functional_position pos_separator functional_position SCB;

fahfa_description : fahfa_pos DASH fahfa DASH | fahfa_pos DASH ROB fahfa RCB DASH;
fahfa : fatty_acid;
fahfa_pos : number;

/* separators */
SPACE : ' ';
COLON : ':';
SEMICOLON : ';';
DASH : '-';
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


