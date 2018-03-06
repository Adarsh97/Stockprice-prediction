(setf *random-state* (make-random-state t))
(setf ActWt (make-array '(5))) (dotimes (j 5) (setf (aref ActWt j) 1) )
(setf test (make-array '(5))) (dotimes (j 5) (setf (aref test j) 0) )
(setf temparray (make-array '(14)))
(setf PERM_VAL .25)
(setf m 0)
(setf n 0)



(setf VARIANCE 1414)

(setf PS_BEST (make-array '(5))) (dotimes (i 5) (setf (aref PS_BEST i) 1) )
(setf AW_BEST (make-array '(5))) (dotimes (i 5) (setf (aref AW_BEST i) 1) )

(defun Get_Variance(ActWt DS PS)

(setf testSize (truncate (* m .3)))

(setf z 0)
(setf y 0)
(setf storage 0)
(setf fstorage 0)
(setf tVARIANCE 0)

(dotimes (z testSize)
(setf fstorage 0)
(dotimes (y 5)
(setf storage (* (aref ActWt y) (expt (aref DS z y) (aref PS y)) ) )
(setf fstorage (+ fstorage storage))
)
(setf tVARIANCE (+ tVARIANCE (expt (- (aref DS z 5) fstorage) 2) ))

)
(return-from Get_Variance tVARIANCE)
)

(defun main()

(defparameter *p* (open "myfile.txt"))
(setf m (read *p*))
(setf n (read *p*))
(setf g (list m n ))
(setf DS (make-array g))
(dotimes (r m)
(dotimes (s n)
(setf (aref DS r s ) (read *p*))
)
)
(close *p*)

(normalize DS)
(setf AW (make-array '(14 5)))
(dotimes (i 14)
   (dotimes (j 5)
      (setf (aref AW i j) (/ (- (rand-val 2415) 1000) 1000))
   )
)

(setf PS (make-array '(5)))
(dotimes (i 5) (setf (aref PS i) 1))
(setf OldActWt (make-array '(5)))
(dotimes (j 5) (setf (aref OldActWt j) (aref ActWt j)))

(dotimes (i 47) ;;47
(dotimes (j 5) (setf (aref OldActWt j) (aref ActWt j)) )
(if (< VARIANCE PERM_VAL) (return))
(setq AW (genetic AW (2d-array-to-list DS) (list (aref PS_BEST 0)(aref PS_BEST 1) (aref PS_BEST 2) (aref PS_BEST 3) (aref PS_BEST 4))))
(Attrbt_Selector AW ActWt)

(setf VARIANCE (Get_Variance ActWt DS PS))
;;(Get_Limiter DS PERM_VAL)
;;(Best_Sofar VARIANCE PERM_VAL)
(Modify_PS PS ActWt OldActWt)
;;(print PS)
(print ActWt)
(print PS)
)

;;(print AW)
;;(print DS)
(checker PS ActWt)

)
(defun checker (PS Actwts)
(defparameter *q* (open "myfile1.txt"))
(setf m (read *q*))
(setf n (read *q*))
(setf g (list m n ))
;;(setf n (- n 1))
(setf AWSET (make-array g))

(dotimes (r m)
(dotimes (s n)
(setf (aref AWSET r s ) (read *q*))
)
)
(close *q*)
(normalize AWSET)
(print AWSET)
(setf sum 0)
(dotimes (i m)
(setf sum 0)
(dotimes (j (- n 1))
(setf sum (+ sum (* (aref Actwt j) (expt (aref AWSET i j) (aref PS_BEST j))))))
(print sum)
(print (aref AWSET i 5))
(terpri)
)
)
(defun Modify_PS(PS ActWt OldActWt)
(dotimes (j 5)
(setf (aref OldActWt j) (- (aref OldActWt j) (aref ActWt j)))
)
(setf max -9999)
(setf min 9999)
(setf minpos 0)
(setf maxpos 0)
(dotimes (j 5)
(setf oldmin min)
(if (< (aref OldActWt j) min) (setf min (aref OldActWt j)))
(if (/= oldmin min) (setf minpos j))
(setf oldmax max)
(if (> (aref OldActWt j) max) (setf max (aref OldActWt j)))
(if (/= oldmax max) (setf maxpos j))
)
(if (/= (aref PS minpos) 1) (setf (aref PS minpos) (- (aref PS minpos) 1)))
(setf (aref PS maxpos) (+ (aref PS maxpos) 1))
)

(defun rand-val(tr)
(+ 0.0 (random tr))
)

(defun normalize (DS)
(setf max 0)
(setf min 999999)
(setf h1 0)
(dotimes (h1 m)
(if (< (aref DS h1 4) min)
(setf min (aref DS h1 4))
)
(if (> (aref DS h1 4) max)
(setf max (aref DS h1 4))
)
)
(setf h1 0)
(dotimes (h1 m)
(setf (aref DS h1 4) (/ (- (aref DS h1 4) min) (- max min) ))
)
)

(defun Attrbt_Selector (AW ActWt)
(setf h1 0)
(setf h2 0)
(dotimes (h1 5)
(dotimes (h2 14)
(setf (aref temparray h2) (aref AW h2 h1))
)
(setf temparray (sort temparray #'<))
(setf medn (/ (+ (aref temparray 6) (aref temparray 7)) 2))
(setf h2 0)
(setf mean 0)
(dotimes (h2 14)
(setf mean (+ mean (aref temparray h2)))
)
(setf mean (/ mean 14))
(setf (aref ActWt h1) (- (* medn 3) (* 2 mean)))
)
)

(defun power(base exp)
 ( if (eq exp 0) 1 ( * base (power base (- exp 1)) ) )
)
;(defun fx (list)
;	(return (nth 1 list)))

;(defvar ith-row-aw
;	(fx AW))

;(defvar ith-row-ds
;	(fx AW))
;Calcualtes the weighted sum = ax + by + cz + .. so on
(defun wt-sum (ith-row-aw ith-row-ds powerset)
(if (not (null ith-row-aw))
  (+ (* (car ith-row-aw) (power (car ith-row-ds) (car powerset))) (wt-sum (cdr ith-row-aw) (cdr ith-row-ds) (cdr powerset)))
  0
)
)


(defvar sum 0)

;Calculates % closeness btw 2 rows of AW matrix
(defun closeness (ith-row-aw dataset powerset)
;(let (sum-value (wt-sum (fx AW) (fx dataset) powerset))
;	(cl-val 0)
  (setq sum 0)
;	((/ (abs (- sum-value (nth 6 (fx dataset)))) (nth 6 (fx dataset))))
;)
(loop for x in dataset
  do (setq sum (+ sum (/ (abs (- (wt-sum ith-row-aw x powerset) (nth 5 x))) (nth 5 x)) ))
  )
  (/ sum 13)
)


;	(loop for x in dataset
;	do (setq sum (+ sum (wt-sum ith-row-aw x)))
;	)
;	sum
;)

(defun 2d-array-to-list (array)
(loop for i below (array-dimension array 0)
      collect (loop for j below (array-dimension array 1)
                    collect (aref array i j))))

(defun list-to-2d-array (list1)
(make-array (list (length list1)
                  (length (first list1)))
            :initial-contents list1))

(defvar fitList '())

(defun fitness (AW dataset powerset)
(setf fitList '())
(loop for x in (2d-array-to-list AW)
  do (
    setq fitList (append fitList (list (/ 1 (+ 1 (closeness x  dataset powerset) ))))
  )
)
;(/ 1 (+ 1 (closeness (list (aref AW i 0) (aref AW i 1) (aref AW i 2) (aref AW i 3) (aref AW i 4))  dataset powerset) ))

)

(setf totalsum 0.0)
(setf partsum 0.0)
(setf count 0)
(setf i 0)
(defun roulette-wheel (fitList)
(setf totalsum (loop for x in fitList sum x))
(setf r (rand-val totalsum))
(setf partsum 0.0)
(setf count 0)
(loop for x in fitList
  do(progn

      (if (> partsum r)
      (setf i count)
      (progn
        (setf partsum (+ partsum x))
        (incf count)
      )
    )

  )
)
count
)

(defun get-list (i aw_list)
(if (= i 1)
  (car aw_list)
  (get-list (decf i) (cdr aw_list))
)
)

(setf newpop '())
(setf p1 '())
(setf p2 '())
(setf c2 '())
(setf c1 '())

(defun cross-over (AW fitList)
(setf newpop '())
(dotimes (i 3) (progn
(setf p1 (get-list (roulette-wheel fitList) (2d-array-to-list AW)))
(setf p2 (get-list (roulette-wheel fitList) (2d-array-to-list AW)))
(setf c1 '())
(setf c2 '())
(loop for y from 1 to 5
  do(progn
    (if (oddp y)
      ;uniform crossover
      (progn
              (setf c1 (append c1 (list (car p1))))
              (setf c2 (append c2 (list (car p2))))
              )


    (progn
      (setf c1 (append c1 (list (car p2))))
      (setf c2 (append c2 (list (car p1))))
    )
    )
    (setf p1 (cdr p1))
    (setf p2 (cdr p2))
    )

  )

  (setf newpop (append newpop (list c1)))
  (setf newpop (append newpop (list c2)))
  )

)
)



(defvar newfitList '())

(defun fitness-again (aw_list dataset powerset)
(setq newfitList '())
(loop for x in aw_list
  do(
    setf newfitList (append newfitList (list (/ 1 (+ 1 (closeness x  dataset powerset) ))))
  )
)
)

(defun maxq (list1)
(if (eq (length list1) 1)
  (car list1)
  (if (> (car list1) (maxq (cdr list1)))
    (car list1)
    (maxq (cdr list1))
  )
)
)

(defun minq (list1)
(if (eq (length list1) 1)
  (car list1)
  (if (< (car list1) (minq (cdr list1)))
    (car list1)
    (minq (cdr list1))
  )
)
)

(defun remove-min (list1)
(
  if (eq (length list1) 1)
  '()
  (
    if (= (car list1) (minq list1))
    (cdr list1)
    (append (list (car list1)) (remove-min (cdr list1)))
  )
)
)

(defun remove-max (list1)
(
  if (eq (length list1) 1)
  '()
  (
    if (= (car list1) (maxq list1))
    (cdr list1)
    (append (list (car list1)) (remove-max (cdr list1)))
  )
)
)

(defun remove-index (i list1)
(
  if(= i 1)
  (cdr list1)
  (append (list (car list1)) (remove-index (decf i) (cdr list1)))
)
)

(defun search-index (val list1)
(
  if(= val (car list1))
  1
  (+ 1 (search-index val (cdr list1)))
)
)
(defun replace-chromes (aw_list newpop fitList newfitList)
      (if  (and (not (null newfitList)) (> (maxq newfitList) (minq fitList)))
				(progn (setq aw_list (remove-index (search-index (minq fitList) fitList) aw_list))
				(setq aw_list (append aw_list (list (get-list (search-index (maxq newfitList) newfitList) newpop))))
				(setq newpop (remove-index (search-index (maxq newfitList) newfitList) newpop))
				(setq fitList (remove-index (search-index (minq fitList) fitList) fitList))
				(setq newfitList (remove-index (search-index (maxq newfitList) newfitList) newfitList)))
)

(if  (and (not (null newfitList)) (> (maxq newfitList) (minq fitList)))
	(progn (setq aw_list (remove-index (search-index (minq fitList) fitList) aw_list))
	(setq aw_list (append aw_list (list (get-list (search-index (maxq newfitList) newfitList) newpop))))
	(setq newpop (remove-index (search-index (maxq newfitList) newfitList) newpop))
	(setq fitList (remove-index (search-index (minq fitList) fitList) fitList))
	(setq newfitList (remove-index (search-index (maxq newfitList) newfitList) newfitList)))
)
(if  (and (not (null newfitList)) (> (maxq newfitList) (minq fitList)))
	(progn (setq aw_list (remove-index (search-index (minq fitList) fitList) aw_list))
	(setq aw_list (append aw_list (list (get-list (search-index (maxq newfitList) newfitList) newpop))))
	(setq newpop (remove-index (search-index (maxq newfitList) newfitList) newpop))
	(setq fitList (remove-index (search-index (minq fitList) fitList) fitList))
	(setq newfitList (remove-index (search-index (maxq newfitList) newfitList) newfitList)))
)
(if  (and (not (null newfitList)) (> (maxq newfitList) (minq fitList)))
	(progn (setq aw_list (remove-index (search-index (minq fitList) fitList) aw_list))
	(setq aw_list (append aw_list (list (get-list (search-index (maxq newfitList) newfitList) newpop))))
	(setq newpop (remove-index (search-index (maxq newfitList) newfitList) newpop))
	(setq fitList (remove-index (search-index (minq fitList) fitList) fitList))
	(setq newfitList (remove-index (search-index (maxq newfitList) newfitList) newfitList)))
)
(if  (and (not (null newfitList)) (> (maxq newfitList) (minq fitList)))
	(progn (setq aw_list (remove-index (search-index (minq fitList) fitList) aw_list))
	(setq aw_list (append aw_list (list (get-list (search-index (maxq newfitList) newfitList) newpop))))
	(setq newpop (remove-index (search-index (maxq newfitList) newfitList) newpop))
	(setq fitList (remove-index (search-index (minq fitList) fitList) fitList))
	(setq newfitList (remove-index (search-index (maxq newfitList) newfitList) newfitList)))
)
(if  (and (not (null newfitList)) (> (maxq newfitList) (minq fitList)))
	(progn (setq aw_list (remove-index (search-index (minq fitList) fitList) aw_list))
	(setq aw_list (append aw_list (list (get-list (search-index (maxq newfitList) newfitList) newpop))))
	(setq newpop (remove-index (search-index (maxq newfitList) newfitList) newpop))
	(setq fitList (remove-index (search-index (minq fitList) fitList) fitList))
	(setq newfitList (remove-index (search-index (maxq newfitList) newfitList) newfitList)))
)



				aw_list
)

;(fitness AW dataset powerset)
;(print fitList)
;(cross-over AW fitList)
;(print newpop)
;(fitness-again newpop dataset powerset)
(defun genetic (AW dataset powerset)

;(setq newAW (list-to-2d-array (replace-chromes (2d-array-to-list AW) newpop fitList newfitList)))

(dotimes (i 2)
  (fitness AW dataset powerset)
  (cross-over AW fitList)
  (setf newfitList '())
  (fitness-again newpop dataset powerset)
  (setq AW (list-to-2d-array (replace-chromes (2d-array-to-list AW) newpop fitList newfitList)))
  ;(fitness AW dataset)
  ;(print fitList)

  )
  AW
  )

(main)
