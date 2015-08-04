; emacs --batch foo.c -l indent -f save-buffer >& log
(require 'cc-mode)
(c-add-style "my-style"
             '("stroustrup"
               (indent-tabs-mode . nil)        ; use spaces rather than tabs
               (c-basic-offset . 4)            ; indent by four spaces
               (tab-width .  4)
               (c-offsets-alist . ((inline-open . 0)  ; custom indentation rules
                                   (brace-list-open . 0)
                                   (innamespace . [0])
                                   (statement-case-open . +)))))


(c-set-style "my-style")        ; use my-style defined above
(auto-fill-mode)
(indent-region (point-min) (point-max) nil)
(untabify (point-min) (point-max) nil)
(delete-trailing-whitespace)
