default : $(PLOTS)

$(PLOTS) : plot.py mode-plot.py
	-python3 plot.py
	-python3 mode-plot.py
    

clean:
	-rm $(PLOTS) __pycache__

.PHONY: clean default
