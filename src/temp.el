(require 'cc-mode)
(require 'cl)
;; regexps stolen and adapted from cc-menus.el

(defun rargs (str formal_targ concretet)
  (replace-regexp-in-string  "typename\\b" ""
   (replace-regexp-in-string (concat formal_targ "\\b") concretet str t))
  )
(defun find-cls-exp (cls)
  (progn
    (goto-char (point-min))
    (re-search-forward
     (concat
      "^"                                  ; beginning of line is required
      "template[ \t\n]*<\\(class\\|typename\\)" ;;1
      "\\([ \t\n]\\|\\\\\n\\)*"                 ;;2
      "\\([^>]+\\)>"                            ;;3
      "\\([ \t\n]\\|\\\\\n\\)*"                 ;;4
      "\\(class\\|struct\\)[ \t]+"              ;;5
      "\\(" cls "\\)"                    ;;6 cls name
      "\\([ \t\n]\\|\\\\\n\\)*;"         ;
      ) (point-max) t) )
  )
(defun find-cls()
  (buffer-substring-no-properties (match-beginning 5)
                                  (match-end 6)
                                  )
  )

(defun find-cls-targ()
  (buffer-substring-no-properties (match-beginning 3)
                                  (match-end 3)
                                  )
  )

;; find method
(defun find-exp(templ-class cls cm)
  (progn
    (goto-char (point-min))
    (re-search-forward
     (concat
      "^\\<"                            ; line MUST start with word char
      "[^()\n]*"                        ; no parentheses before
      "template[ \t\n]*<\\(class\\|typename\\)" ;;1
      "\\([ \t\n]\\|\\\\\n\\)*"                 ;;2
      "\\([^>]+\\)>"                            ;;3 template placeholder
      "\\([ \t\n]\\|\\\\\n\\)*"                 ;;4
      "\\(static\\)?"                           ;;5
      "\\([ \t\n]\\|\\\\\n\\)*"                 ;;6
      "\\(" "[" c-alpha "]" "[" c-alnum " _:<>]+\\)"                ;7 (return type) match identifier chars
      "\\([ \t\n]\\|\\\\\n\\)*"         ;8
      cls
      (if templ-class "<\\([^>]+\\)>" "") ;9
      "::"
      cm                                ; match function name
      "\\([ \t\n]\\|\\\\\n\\)*("        ;9,10 see above, BUT the arg list
      "\\([ \t\n]\\|\\\\\n\\)*"         ;10,11 must not start
      "\\([^ \t\n(*]"                   ;11,12 with an asterisk or parentheses
      "[^()]*\\(([^()]*)[^()]*\\)*"    ; Maybe function pointer arguments
      "\\)?)"
      "\\([ \t\n]\\|\\\\\n\\)*[^ \t\n;(]"
      ) (point-max) nil)
    )
  )
;; find return type
(defun find-ret()
  (buffer-substring-no-properties (match-beginning 7)
                                  (match-end 7)
                                  )
  )

(defun find-args(templ-class)
  (if templ-class
      (buffer-substring-no-properties (match-beginning 12)
                                      (match-end 12)
                                      )
    (buffer-substring-no-properties (match-beginning 11)
                                    (match-end 11)
                                    ))
  )
(defun find-cls-targ-in-method()
  (buffer-substring-no-properties (match-beginning 9)
                                  (match-end 9)
                                  )
  )

(defun find-targ()
  (buffer-substring-no-properties (match-beginning 3)
                                  (match-end 3)
                                  )
  )

(defun instantiate-templates(cls ts-cls cms ts-methods)
  """instantiate a template class and its methods"""
  (interactive)
  ;; clear the inst file
  (with-temp-file (concat cls "_inst.cpp")
    (insert (concat "//-*-mode:c++-*-\n" "//Generated on: " (current-time-string) 
                    ". Do not edit\n")
            )
    )
  ;; now instantiate each method that is passed for each of the template
  ;; args in ts
  (setq templ-class t)
  (save-excursion
    (if (find-cls-exp cls)
        (instantiate-template/class cls (find-cls) (find-cls-targ) ts-cls)
      (setq templ-class nil)
      )
    (dolist (cm cms)
      (progn
        (if (find-exp templ-class cls cm)
            (if templ-class
                (instantiate-template-template/method cls ts-cls cm (find-ret) 
                                                      (find-args templ-class) (find-targ) 
                                                      (find-cls-targ-in-method) ts-methods)
              (instantiate-template/method cls cm (find-ret)
                                           (find-args templ-class) (find-targ) ts-methods)
              )
          (message (format "failed to find %s::%s " cls cm))
          )
        )
      )
    )
  )

(defun instantiate-template/class(cls-small cls-expanded formal_targ targs)
  (with-temp-buffer
    (message "instantiating cls: " cls-small)
    (insert (concat "//________________________________________________________\n"
                    "// class:" cls-small "\n"))
    (dolist (targ targs)
      (insert (format "template %s<%s>;\n"  cls-expanded targ ))
      )
    (append-to-file (point-min) (point-max) (concat cls-small "_inst.cpp") )
    )
  )

(defun instantiate-template/method(cls curmethod ret args formal_targ targs)
  """instantiate a template function"""
  (let ((curcls cls)
        (fname cls)
        )
    (with-temp-buffer
      (message "instantiating: " cm)
      (insert (concat "//________________________________________________________\n"
                      "// method:" curmethod "\n"))
      (dolist (targ targs)
        (insert (format "template %s %s::%s(%s);\n" (rargs ret formal_targ targ)
                        curcls cm (rargs args formal_targ targ ) ) )
        )
      (append-to-file (point-min) (point-max) (concat fname "_inst.cpp") )
      )
    )
  )

(defun instantiate-template-template/method(cls ctargs curmethod ret args formal_targ cls-formal-targ targs)
  """instantiate a template function"""
  (let ((curcls cls)
        (fname cls)
        )
    (with-temp-buffer
      (message "instantiating: " cm)
      (insert (concat "//________________________________________________________\n"
                      "// method:" curmethod "," cls-formal-targ "," formal_targ "\n"))
      (dolist (ctarg ctargs)
        (dolist (targ targs)
          (insert (format "template %s %s<%s>::%s(%s);\n" 
                          (rargs (rargs ret formal_targ targ) cls-formal-targ ctarg)
                          curcls ctarg cm
                          (rargs (rargs args formal_targ targ ) cls-formal-targ ctarg)
                          ) )
          ))
      (append-to-file (point-min) (point-max) (concat fname "_inst.cpp") )
      )
    )
  )
