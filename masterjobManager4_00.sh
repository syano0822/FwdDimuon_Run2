#!/bin/sh

#############################################################################
# Ver 4.00 (2015/11/22) by terakoya.k                                       #
#                                                                           #
# To manage all subjobs at the same time                                    #
# Source this file in aliensh after selecting option (in "Select mode")     #
#############################################################################


# Select mode (Please select mode by yourself) ------------------------------

# Display color mode (1: color, other: monochrome)
colorMode=1

# Selection mode
# (This script manages... 1: masterjob(s) you select, other: all jobs you own)
selectMode=0

# Yes/No mode (1: confirm action before doing it, other: not confirm)
# [1/other] is recommended for [selection mode/all job mode] (i.e. selectMode=confirmMode is recommended)
confirmMode=0

# ---------------------------------------------------------------------------


# Color mode ----------------------------------------------------------------
if [ ${colorMode} -eq 1 ]; then
	Cm="\033[35m" # Magenta
	Cr="\033[31m" # Red
	Cy="\033[33m" # Yellow
	Cb="\033[34m" # Blue
	Cc="\033[36m" # Cyan
	Cz="\033[m"   # Normal
else
	Cm="\033[m"
	Cr="\033[m"
	Cy="\033[m"
	Cb="\033[m"
	Cc="\033[m"
	Cz="\033[m"
fi
# ---------------------------------------------------------------------------


# Check if you are in aliensh -----------------------------------------------
# by checking if "masterjob" command exists or not (Is there better way?)
typeOutput=`type -t masterjob` # OUTPUT: "alias" -> OK, "" -> bad

if [[ ${typeOutput} = "alias" ]]; then
	# you are in aliensh
	isMain1=1 # go to next step
	echo
else
	# you are not in aliensh
	isMain1=0 # break
	echo
	echo -e "${Cm}   ### You are not in aliensh. Source this file after doing \"aliensh\".${Cz}"
	echo
fi
# ---------------------------------------------------------------------------


# Display option ------------------------------------------------------------
if [[ ${isMain1} == 1 ]]; then
	echo
	echo -e "${Cb}   # Select number from below. (Version 4.00)${Cz}"
	echo -e "${Cb}   (ATTENTION! Number changed since ver3.00.)${Cz}"
	echo
	echo -e "${Cb}   +------+---------------------------------+${Cz}"
	echo -e "${Cb}   |  No. |             ACTION              |${Cz}"
	echo -e "${Cb}   +------+---------------------------------+${Cz}"
	echo -e "${Cb}   |   ${Cy}0. ${Cb}|  ${Cy}SHOW     all        MASTERjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |   ${Cr}1. ${Cb}|  ${Cr}KILL     all DONE      SUBjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |   ${Cc}2. ${Cb}|  ${Cc}RESUBMIT all ERROR     SUBjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |   ${Cc}3. ${Cb}|  ${Cc}RESUBMIT all ZOMBIE    SUBjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |   ${Cc}4. ${Cb}|  ${Cc}RESUBMIT all E_SPLT MASTERjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |   ${Cr}5. ${Cb}|  ${Cr}KILL     all DONE   MASTERjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |   ${Cm}7. ${Cb}|  ${Cm}KILL     all ERROR     SUBjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |   ${Cm}8. ${Cb}|  ${Cm}KILL     all ZOMBIE    SUBjobs ${Cb}|${Cz}"
	echo -e "${Cb}   |  ${Cm}99. ${Cb}|  ${Cm}KILL     ALL        MASTERjobs ${Cb}|${Cz}"
	echo -e "${Cb}   +------+---------------------------------+"
	if [[ ${selectMode} == 1 ]]; then
		echo -e "${Cb}   Selection Mode:${Cr} ON${Cb} (you select job numbers later)${Cz}"
		echo -e "${Cb}   (You can change this option in the script)${Cz}"
	else
		echo -e "${Cb}   Selection Mode:${Cr} OFF${Cb} (manage all jobs you own)${Cz}"
		echo -e "${Cb}   (You can change this option in the script)${Cz}"
	fi
	echo
# ---------------------------------------------------------------------------


# Input option number -------------------------------------------------------
	echo -ne "${Cb}Your choice (type a number and press ENTER): ${Cz}"
	read nMode

	echo
# ---------------------------------------------------------------------------


# Show mode you selected ----------------------------------------------------
	isMain2=1 # default
	case "${nMode}" in
		0) echo -e "${Cy}   0: SHOW MASTERjobs${Cz}"
				echo
				ps -mA
				isMain2=0;; # break
		1)  echo -e "${Cr}   1:  \"KILL\" all \"DONE\" SUBjobs${Cz}";;
		2)  echo -e "${Cc}   2:  \"RESUBMIT\" all \"ERROR\" SUBjobs${Cz}";;
		3)  echo -e "${Cc}   3:  \"RESUBMIT\" all \"ZOMBIE\" SUBjobs${Cz}";;
		4)  echo -e "${Cc}   4:  \"RESUBMIT\" all \"ERROR_SPLIT\" MASTERjobs${Cz}";;
		5)  echo -e "${Cr}   5:  \"KILL\" all \"DONE\" MASTERjobs${Cz}";;
		7)  echo -e "${Cm}   7:  \"KILL\" all \"ERROR\" SUBjobs${Cz}";;
		8)  echo -e "${Cm}   8:  \"KILL\" all \"ZOMBIE\" SUBjobs${Cz}";;
		99) echo -e "${Cm}   99: \"KILL\" \"ALL\" MASTERjobs${Cz}";;
		*)  echo -e "${Cb}   ### Nothing to do${Cz}"
				isMain2=0;; # barek
	esac
	echo
# ---------------------------------------------------------------------------


# Select PID (masterjob number) ---------------------------------------------
	if [[ ${isMain2} == 1 && ${selectMode} == 1 ]]; then
		while [[ 1 ]]; do
			echo -ne "${Cb}Management starts from (type a 9-digit number and press ENTER): ${Cz}"
			read nStart
			if [[ ${nStart} =~ ^[0-9]{9}$ ]]; then
				break
			else
				echo -e "${Cm}   ### Your input is not 9-digit number.${Cz}"
			fi
		done

		while [[ 1 ]]; do
			echo -ne "${Cb}Management ends in (type a 9-digit number and press ENTER):     ${Cz}"
			read nEnd
			if [[ ${nEnd} =~ ^[0-9]{9}$ ]]; then
				break
			else
				echo -e "${Cm}   ### Your input is not 9-digit number.${Cz}"
			fi
		done

	fi

# ---------------------------------------------------------------------------


# Main ----------------------------------------------------------------------
	if [[ ${isMain2} == 1 ]]; then
		isMain3=1 # default
		if [[ ${confirmMode} == 1 ]]; then
			if [[ ${selectMode} == 1 ]]; then
				echo
				echo -ne "${Cb}   Is it OK? (From: ${Cz}"
				echo -n ${nStart}
				echo -ne "${Cb} To: ${Cz}"
				echo -n ${nEnd}
				echo -e "${Cb})${Cz}"
			else
				echo
				echo -e "${Cb}   Is it OK? (Managing all jobs you own)${Cz}"
			fi

			while [[ 1 ]]; do
				echo -ne "${Cb}   (y/n): ${Cz}"
				read yesno
				if [[ ${yesno} = "n" || ${yesno} = "no" ]]; then
					isMain3=0 # break
					echo
					echo -e "${Cb}   ### Nothing to do${Cz}"
					echo
					break
				elif [[ ${yesno} = "y" || ${yesno} = "yes" ]]; then
					echo
					break
				else
					echo -e "${Cm}   ### Please type "y"(Yes) or "n"(No).${Cz}"
				fi
			done
		fi # if [[ ${confirmMode} == 1 ]]; then

		if [[ ${isMain3} == 1 ]]; then
			ps -mA | if [[ ${nMode} == 4 ]]; then awk '$5~"ESP"{print $2}'; elif [[ ${nMode} == 5 ]]; then awk '$6~"D"{print $2}'; else awk '{print $2}'; fi | /usr/bin/col -bx | sed -e "s/.\{1,5\}m//g" | while read line

			# "ps -mA"     : show the jobidlist(-m) which you own(-A)
			# "awk ..."    : remove username and get jobidlist from ps output
			# "col ..."    : remove color codes 1 (e.g. ESC[0mESC[1m) from awk output
			# "sed -e ..." : remove color codes 2 (e.g. ESC[0mESC[1m) from awk output
			do
				if [[ ${selectMode} == 1 ]]; then
					if [[ ${line} < ${nStart} || ${line} > ${nEnd} ]]; then
						continue
					fi
				fi
				echo
				echo -ne "For masterjob: "
				echo -e ${line}
				case "${nMode}" in
					1) masterjob ${line} -status DONE      kill;;
					2) masterjob ${line} -status ERROR_ALL resubmit;;
					3) masterjob ${line} -status ZOMBIE resubmit;;
					4) resubmit  ${line};;
					5) kill      ${line};;
					7) masterjob ${line} -status ERROR_ALL kill;;
					8) masterjob ${line} -status ZOMBIE kill;;
					99) kill     ${line};;
				esac

			done

		fi # if [ ${isMain3} == 1 ]; then

	fi # if [ ${isMain2} == 1 ]; then

fi # if [ ${isMain1} == 1 ]; then
# ---------------------------------------------------------------------------
